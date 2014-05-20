#!/usr/bin/env python

from Cantera import *
import numpy as np

class flow():

    def __init__(self):

        self.gas = GRI30()
        self.gas.setMoleFractions('N2:0.79,O2:0.21')

        self.R = self.gas.cp_mass() - self.gas.cv_mass()
        self.gam = self.gas.cp_mass()/self.gas.cv_mass()

        self.a = np.sqrt(self.gam*self.R*self.gas.temperature())

        self.M = 2.0

        self.v = self.M*self.a
        
        self.setMassFlux()
        self.setMomentumFlux()
        self.setTotalEnthalpy()

    def setTotalEnthalpy(self):

        self.Ht = self.gas.enthalpy_mass() + self.v**2/2.0

    def setMomentumFlux(self):

        self.momFlux = self.gas.pressure() + self.gas.density()*self.v**2

    def setMassFlux(self):

        self.massFlux = self.gas.density()*self.v
    
