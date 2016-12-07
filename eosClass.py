# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 14:52:51 2015

@author: Vishwanath Haily Dalvi
"""
import scipy, scipy.optimize, scipy.integrate

def funcPsat(P, T, molecule, alpha):
    if abs(P.imag) > 1e-9:
        raise RuntimeError('Imaginary P found')
    else:
        P = P.real
    [ZL, ZG] = molecule.Z(T, P, alpha = alpha)
        
    gRL = molecule.gR(T, P, ZL, alpha = alpha)
    gRG = molecule.gR(T, P, ZG, alpha = alpha)
    return (gRG-gRL)/(molecule.R*T)


class EOS:
    """
    Fives types of EOS supported: srk, pr, pt, hkm2
    """
    def __init__(self, molecule):
        '''
        mol = EOS(molecule)
        Here: molecule is of class "Compound" from Perry's data
        '''
        self.R = 8314.0 #J/kmol-K
        self.setMolecule(molecule)   
    def setMolecule(self, molecule):
        self.Molecule = molecule
        self.Name = molecule.Name
        self.Tc = molecule.Tc
        self.Pc = molecule.Pc
        self.acc = molecule.acc
        self.Zc = molecule.Zc   
    def setEOS(self, typeeos, OmegaB = False):
        self.a, self.b, self.c, self.d = self.getabcd(typeeos, OmegaB = OmegaB)
        self.typeeos = typeeos
    def getOmegas(self, typeeos, OmegaB = False):
        acc = self.acc
        if typeeos in ['srk']:
            OmegaC = 0.0
            onethird = 1.0/3.0
            OmegaB = onethird*(2**onethird-1)
            OmegaD = OmegaB + 0.0
            OmegaA = onethird**2/(2**onethird-1)
        elif typeeos == 'pr':
            OmegaB = 0.077796 
            OmegaA = 0.457235
            OmegaC = OmegaB + 0.0
            OmegaD = OmegaB + 0.0
        elif typeeos == 'pt':
            Zc = 0.3272 - 0.0537*acc - 0.0147*acc**2 #Note: NOT experimental Zc!
            A = 1.0
            B = 2 - 3*Zc
            C = 3*Zc**2
            D = -Zc**3
            rts = scipy.roots([A, B, C, D])
            OmegaB = min([r.real for r in rts if (r.real > 0.0) and (abs(r.imag) < 1e-9)])
            OmegaD = OmegaB + 0.0
            OmegaC = 1 - 3*Zc
            OmegaA = 3*Zc**2 + OmegaB*OmegaC + OmegaC*OmegaD + OmegaB*OmegaD + OmegaC + OmegaD
        elif typeeos in  ['hkm2']:
            Zc = 0.3175 - 0.0364*acc -0.0245*acc**2  #Note: NOT experimental Zc!
            A = 1.0
            B = 1.25 - 3*Zc
            C = 3*Zc**2 - 1.5*Zc + 0.5
            D = -Zc**3
            rts = scipy.roots([A, B, C, D])
            OmegaB = min([r.real for r in rts if (r.real > 0.0) and (abs(r.imag) < 1e-9)])
            OmegaE = 6*Zc - 3*OmegaB - 2
            OmegaA = 3*Zc**2 - 0.5*OmegaB**2 - 0.75*OmegaB*OmegaE - 0.5*(OmegaB+OmegaE)
            n = -0.5; m = -0.5
            nbplusme = n*OmegaB + m*OmegaE
            nbme = n*OmegaB*m*OmegaE
            OmegaD = 0.5*(nbplusme - scipy.sqrt(nbplusme**2+4*nbme))
            OmegaC = -nbme/OmegaD
        return OmegaA, OmegaB, OmegaC, OmegaD
    def getabcd(self, typeeos, OmegaB = False):
        RTcPc = self.R*self.Tc/self.Pc
        RTc2Pc = RTcPc*self.R*self.Tc
        OmegaA, OmegaB, OmegaC, OmegaD = self.getOmegas(typeeos, OmegaB = OmegaB)
        a = RTc2Pc*OmegaA.real
        b = RTcPc*OmegaB
        c = RTcPc*OmegaC
        d = RTcPc*OmegaD
        return a, b, c, d
    def alpha(self, T):
        Tr = T/self.Tc
        typeeos = self.typeeos
        acc = self.acc
        if typeeos == 'pr':
            k = 0.37464 + 1.54226*acc - 0.26992*acc**2
            sqrtalpha = 1 + k*(1-Tr**0.5)
            alpha = sqrtalpha**2
        elif typeeos == 'srk':
            k = 0.48 + 1.574*acc - 0.176*acc**2
            sqrtalpha = 1 + k*(1-Tr**0.5)
            alpha = sqrtalpha**2
        elif typeeos == 'pt':
            F = 0.452413+ 1.30982*acc -0.295937*acc**2
            alpha = (1 + F*(1-Tr**0.5))**2
        elif typeeos == 'hkm2':
            k1 = 0.0821 + 0.3042*acc - 0.0730*acc**2
            k2 = 3.058 + 1.5479*Tr
            alpha = scipy.exp(k2*(1-Tr**k1))
        else:
            raise ValueError("Alpha function for %s not yet defined"%typeeos)
        return alpha

    def P(self, T, v, alpha=False):
        if not alpha:
            alpha = self.alpha(T)
        a, b, c, d = self.a*alpha, self.b, self.c, self.d
        P = self.R*T/(v-b) - a/(v**2 + (c+d)*v - c*d)
        return P.real
    def Z(self, T, P, alpha = False):
        P = 1e-20 if P == 0.0 else P
        if not alpha:
            alpha = self.alpha(T)
        PRT = P/(self.R*T)
        a, b, c, d = self.a*alpha, self.b, self.c, self.d
        A = PRT*(c + d - b) - 1
        B = PRT**2*(a/P - b*c - c*d - d*b) - PRT*(c+d)
        C = PRT**3*(b*c*d+c*d/PRT-a*b/P)
        roots = scipy.roots([1, A.real, B.real, C.real])
        roots = [x.real for x in roots if abs(x.imag) < 1e-16]
        if len(roots) == 3:
            return [min(roots), max(roots)]
        elif len(roots) == 1:
            return [roots[0]]
        else:
            return []            
    def gR(self, T, P, Z, alpha = False):
        P = 1e-20 if P == 0.0 else P
        if not alpha:
            alpha = self.alpha(T)
        a, b, c, d = self.a*alpha, self.b, self.c, self.d
        disc = (c+d)**2 + 4*c*d; disc = disc.real
        A = 0.5*(c+d+scipy.sqrt(disc))
        B = 0.5*(c+d-scipy.sqrt(disc))
        PRT = P/(self.R*T)
        gRbyRT = Z - 1.0 - scipy.log(Z-b*PRT) - a*PRT/P/(A-B)*scipy.log((Z+A*PRT)/(Z+B*PRT))
        return self.R*T*gRbyRT.real
    def getvborder(self, T, alpha = False):
        if not alpha:
            alpha = self.alpha(T)
        RT = self.R*T
        a, b, c, d = self.a*alpha, self.b, self.c, self.d
        a4 = -1.0
        a3 = -2*(c+d) + 2*a/RT
        a2 = -c**2 - d**2 - 4*a*b/RT + a*c/RT + a*d/RT
        a1 = 2*c**2*d + 2*c*d**2 + 2*a*b**2/RT - 2*a*b*c/RT - 2*a*b*d/RT
        a0 = -c**2*d**2 + a*b**2*c/RT + a*b**2*d/RT
        roots = scipy.roots([a4.real, a3.real, a2.real, a1.real, a0.real])
        roots = [x.real for x in roots if abs(x.imag) < 1e-16 and x.real > b]
        if len(roots) == 2:
            vborder = [roots[0], roots[1]]
            return True, min(vborder), max(vborder)
        else:
            return False, None, None
    def getPminPmax(self, T, alpha = False):
        if not alpha:
            alpha = self.alpha(T)

        boolean, vmin, vmax = self.getvborder(T, alpha = alpha)
        if boolean:
            Pmin = self.P(T, vmin, alpha = alpha)
            Pmax = self.P(T, vmax, alpha = alpha)
            return True, Pmin, Pmax
        else:
            return False, None, None
            
    def Psat(self, T, alpha = False):
        if not alpha:
            alpha = self.alpha(T)
        Tc = self.Tc
        if T > Tc:
            return False, None
        elif T == Tc:
            return True, self.Pc
        else:
            boolean, Pmin, Pmax = self.getPminPmax(T, alpha = alpha)
            if boolean:
                if Pmin < 0.0: 
                    Pguess = 0.0
                else:
                    Pguess = 0.999*Pmin + 0.001*Pmax
                solution = scipy.optimize.newton(funcPsat, Pguess, args = (T, self, alpha))
                return boolean, solution
            else:
                return boolean, None


        
        
    
    
        
