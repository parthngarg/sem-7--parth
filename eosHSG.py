
# coding: utf-8

# # Inheritance

# It is sometimes nice to have an object that builds on another object.  For example, in our perrysdata.Compound class, we have not programmed for getting ideal gas free energy, enthalpy or entropy.  But these are so commonly used that it would be a good idea to put them in.  Pretend that, for whatever reason, we cannot modify the original class.  So we make a new class into which we <i>inherit</i> the old class.  

# In[4]:

# import scipy and perrydata.Compound
import scipy
import scipy.integrate as scint
from perrysdata import Compound


# In[15]:

#Now lets inherit
class CompoundH(Compound):  #See how the class line takes an argument?  The argument is the parent.  This class is the child. 
    def __init__(self, name):
        Compound.__init__(self, name)  #This calls the constructor of the Parent class. 
    def Hig(self, T):  #Enthalpy of ideal gases is a function of T alone (this is not in Parent class)
        Tf = self.Tf #Reference temperature (from Parent Class)
        Hf = self.Hf #Enthalpy of formation (from Parent Class)
        Cp = self.CpIG #Ideal gas heat capacity function from Parent class (not how we assigned the function to a pointer)
        Hig = Hf + scint.quad(Cp, Tf, T)[0] #Integrate Cp from Tf to T.  Note: Quad returns a 2 element tuple.  The first element is the integral
        return Hig


# In[29]:

#Lets see if it works!  
if __name__ == "__main__": #Guess why this is here!
    methane = CompoundH("Methane")
    Tf = methane.Tf  #Old parent class
    Hf = methane.Hf
    Cp = methane.CpIG
    T = 180.0 #K
    Hig = Hf + scint.quad(Cp, Tf, T)[0]
    print 'Hig = {} J/kmol'.format(Hig)
    print 'Hig from class = {} J/kmol'.format(methane.Hig(T))


# In[22]:

#Great, lets progam in the function for Sig and Gig.  But lets not modify the above class.  Lets inherit some more.
class CompoundHSG(CompoundH): #Awwwww! Our little child is all grown up!!  It has a child of it own. How cute!!!
    def __init__(self, name):
        CompoundH.__init__(self, name)
    def Sig(self, T, P):  #Entropy needs a pressure!
        Tf = self.Tf
        Pf = self.Pf
        Sf = self.Sf
        Cp = self.CpIG
        R = 8314.0 #J/kmol-K
        Sig = Sf + scint.quad(lambda T:Cp(T)/T, Tf, T)[0] - R*scipy.log(P/Pf)
        return Sig
    def Gig(self, T, P):
        Hig = self.Hig(T) #our inheritance!
        Sig = self.Sig(T, P)
        Gig = Hig - T*Sig
        return Gig


# In[30]:

if __name__ == "__main__":
    methane = CompoundHSG("Methane")
    T = 180.0#K
    P = 1.013e5 #Pa
    Hig = methane.Hig(T)
    Sig = methane.Sig(T, P)
    Gig = methane.Gig(T, P)
    print "Ideal gas Enthalpy, Entropy and Free Energy as {} MJ/kmol, {} kJ/kmol and {} MJ/kmol respectively.".format(Hig/1e6,
                                                                                                                      Sig/1e3,
                                                                                                                      Gig/1e6)
    print " "
    print "Ideal gas Enthalpy, Entropy and Free Energy as %.2f MJ/kmol, %.2f kJ/kmol and %6.2f MJ/kmol respectively."%(Hig/1e6,
                                                                                                                      Sig/1e3,
                                                                                                                      Gig/1e6)


# ## Modifying the EOS class
# Incorporating this new CompoundHSG class in the EOS class is simplicity itself.

# In[28]:

from eosClass import EOS as EOS_Parent


# In[35]:

if __name__ == "__main__":
    methane = CompoundHSG("Methane")
    methane = EOS_Parent(methane) #See!  That easy!
    T = 180.0#K
    P = 1.013e5 #Pa

    Hig = methane.Molecule.Hig(T)    #Here we see an attribution chain
    Sig = methane.Molecule.Sig(T, P)
    Gig = methane.Molecule.Gig(T, P)
    print "Ideal gas Enthalpy, Entropy and Free Energy as %.2f MJ/kmol, %.2f kJ/kmol and %6.2f MJ/kmol respectively."%(Hig/1e6,
                                                                                                                      Sig/1e3,
                                                                                                                      Gig/1e6)
    print " "
    methane.setEOS('srk')
    print "Vapour pressure from correlation is %.2f bar but from srk-EOS is %.2f bar"%(methane.Molecule.Pvap(T)*1e-5,
                                                                                     methane.Psat(T)[1]*1e-5)


# But lets <i>modify</i> the EOS class itself to calculate non-ideal enthalpies, entropies and free energies

# In[40]:

import scipy.misc as scmisc
class EOS(EOS_Parent):
    def __init__(self, molecule):
        EOS_Parent.__init__(self, molecule)
    def sR(self, T, P, Z):
        P = 1e-20 if P == 0.0 else P    
        dalphabydT = scmisc.derivative(self.alpha, T, dx = 1e-6) 
        a, b, c, d = self.a, self.b, self.c, self.d
        disc = (c+d)**2 + 4*c*d; disc = disc.real
        A = 0.5*(c+d+scipy.sqrt(disc))
        B = 0.5*(c+d-scipy.sqrt(disc))
        R = self.R
        PRT = P/(R*T)
        S = R*scipy.log(Z - b*PRT) + a*dalphabydT/(A - B)*scipy.log((Z+A*PRT)/(Z+B*PRT))
        return S.real
    def hR(self, T, P, Z):
        P = 1e-20 if P == 0.0 else P 
        alpha = self.alpha(T)
        dalphabydT = scmisc.derivative(self.alpha, T, dx = 1e-6) 
        a, b, c, d = self.a, self.b, self.c, self.d
        disc = (c+d)**2 + 4*c*d; disc = disc.real
        A = 0.5*(c+d+scipy.sqrt(disc))
        B = 0.5*(c+d-scipy.sqrt(disc))
        R = self.R
        PRT = P/(R*T)
        H = R*T*(Z - 1) - a*(alpha - T*dalphabydT)/(A - B)*scipy.log((Z+A*PRT)/(Z+B*PRT))
        return H.real
    def H(self, T, P, Z):
        hR = self.hR(T, P, Z)
        hig = self.Molecule.Hig(T)
        return hig+hR
    def S(self, T, P, Z):
        sR = self.sR(T, P, Z)
        sig = self.Molecule.Sig(T, P)
        return sig+sR
    def G(self, T, P, Z):
        gR = self.gR(T, P, Z)
        gig = self.Molecule.Gig(T, P)
        return gig+gR


# In[50]:

if __name__ == "__main__":
    methane = CompoundHSG("Methane")
    methane = EOS(methane) 
    methane.setEOS('srk')
    T = 180.0#K
    P = methane.Psat(T)[1]
    [ZL, ZG] = methane.Z(T, P)
    print "Residual Liquid Enthalpy, Entropy and Free Energy is %.2f MJ/kmol, %.2f kJ/kmol and %.2f MJ/kmol respectively"%(methane.hR(T, P, ZL)*1e-6,
                                                                                                                           methane.sR(T, P, ZL)*1e-3,
                                                                                                                           methane.gR(T, P, ZL)*1e-6)
    print ""
    print "Residual Vapour Enthalpy, Entropy and Free Energy is %.2f MJ/kmol, %.2f kJ/kmol and %.2f MJ/kmol respectively"%(methane.hR(T, P, ZG)*1e-6,
                                                                                                                           methane.sR(T, P, ZG)*1e-3,
                                                                                                                           methane.gR(T, P, ZG)*1e-6)
    print ""
    Hvap = methane.Molecule.Hvap(T)
    HL = methane.H(T, P, ZL); HG = methane.H(T, P, ZG)
    Hsat = HG - HL
    print "Enthalpy of vapourization by correlation is %.2f MJ/kmol but from %s_EOS is %.2f MJ/kmol"%(Hvap/1e6, 
                                                                                                     methane.typeeos,
                                                                                                     Hsat/1e6)


# # Conclusion
# Now that we are happy with our code, lets write a small "wrapper" to present a useful object to other users.

# In[53]:

def getCompound(name, typeeos):
    comp = CompoundHSG(name)
    comp = EOS(comp)
    comp.setEOS(typeeos)
    return comp


# Now download this as a ".py" file (see Download As in File menu) and place it in this folder.

# In[ ]:



