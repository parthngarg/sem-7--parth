# -*- coding: utf-8 -*-
"""
Created on Tue Jan 20 17:07:29 2015

@author: Ashwin
"""

'''
Perry's Handbook VII:    
Vapour Pressure Table 2-6
Densities of pure substance Table 2-30
Critical Constants and Acentric Factors of Inorganic and Organic Compounds Table 2-164    
Heats of Vapourization:  Table 2-193 
Heat Capacities in Ideal Gas State Table 2-198
Heats and Free Energies of Formation and COmbustion Table 2-221
'''

import os
import pandas
import scipy as sc

filename = 'Data.xlsx'
thisdir = os.getcwd()
if os.path.isdir("Data"):
    os.chdir("Data")
xl_file = pandas.ExcelFile(filename)
os.chdir(thisdir)
df = {}
df["Critical"] = xl_file.parse("Critical")
df["LiqDens"] = xl_file.parse("LiqDens")
df["Pvap"] = xl_file.parse("Pvap")
df["Hvap"] = xl_file.parse("Hvap")
df["CpIG"] = xl_file.parse("CpIG")
df["HandGf"] = xl_file.parse("HandGf")



def getCwater(T):
    if T < 333.15:
        C1, C2, C3, C4 = 5.459, 0.30542, 647.13, 0.081
    elif T >= 333.15 and T < 403.15:
        C1, C2, C3, C4 = 4.9669, 0.27788, 647.13, 0.1874
    else:
        C1, C2, C3, C4 = 4.391, 0.2487, 647.13, 0.2534
    return C1, C2, C3, C4

getCwater = sc.vectorize(getCwater)

class Compound:
    def __init__(self, name):
        self.Name = name
        self.getCriticalConstants()
        self.getHandGofFormation()
        self.constCp = self.getConstants('CpIG')
        self.constPvap = self.getConstants('Pvap')
        self.constHvap = self.getConstants('Hvap')
        self.constLiqDens = self.getConstants('LiqDens')
    def getCriticalConstants(self):
        name = self.Name
        db = df['Critical']
        comp = db[db.Name == name]
        self.MW = comp.MolWt.values[0]  
        self.Tc = comp.Tc.values[0]
        self.Pc = comp.Pc.values[0]
        self.Vc = comp.Vc.values[0]
        self.Zc = comp.Zc.values[0]
        self.acc = comp.Acc.values[0]
    def getHandGofFormation(self):
        name = self.Name
        db = df['HandGf']
        comp = db[db.Name == name]
        self.Hf = comp.Hf.values[0]
        self.Gf = comp.Gf.values[0]
        self.Tf = 298.15 #K
        self.Pf = 1.0e5 #Pa
        self.Sf = (self.Hf - self.Gf)/self.Tf
        self.Sfabinit = comp.Sf.values[0]
        self.Hcomb = comp.Hcomb.values[0]
    def getConstants(self, sheet):
        name = self.Name
        cnsts = ['C1','C2','C3','C4','C5']
        constants = []
        db = df[sheet]
        comp = db[db.Name == name]
        for cnst in cnsts:
            try:
                constants.append(comp[cnst].values[0])
            except KeyError:
                pass
        return constants
    def CpIG(self, T):
        [C1, C2, C3, C4, C5] = self.constCp
        cp  = C1 
        cp += C2*((C3/T)/sc.sinh(C3/T))**2
        cp += C4*((C5/T)/sc.cosh(C5/T))**2
        return cp
    def Pvap(self, T):
        [C1, C2, C3, C4, C5] = self.constPvap
        p = C1 + C2/T + C3*sc.log(T) + C4*T**C5
        return sc.exp(p)
    def Hvap(self, T):
        [C1, C2, C3, C4] = self.constHvap
        Tr = T/self.Tc
        exponent = C2 + C3*Tr + C4*Tr**2
        hvap = C1*(1-Tr)**exponent
        return hvap
    def LiqDens(self, T):
        if self.Name == 'Water':
            C1, C2, C3, C4 = getCwater(T)
        else:
            [C1, C2, C3, C4] = self.constLiqDens    
        exponent = 1 + abs(1-T/C3)**C4
        liqdens = C1/C2**exponent
        return liqdens

def getallnames():
    db = df["Critical"]
    listNames = db.Name.to_dict().values()
    return listNames
    
if __name__ == "__main__":
    met = Compound('Methane')
      
        
