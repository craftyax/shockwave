#!/usr/bin/env python

from Cantera import *
import numpy as np

class flow():

    def __init__(self):

        # Set gas properties
        self.gas = GRI30()
        self.gas.setMoleFractions('N2:0.79,O2:0.21')
        self.R = self.gas.cp_mass() - self.gas.cv_mass()
        self.gam = self.gas.cp_mass()/self.gas.cv_mass()

        # Set flow speed of sound, Mach number, and velocity
        self.a = np.sqrt(self.gam*self.R*self.gas.temperature())
        self.M = 2.0
        self.v = self.M*self.a
        
        # Set conservation variables
        self.setMassFlux()
        self.setMomentumFlux()
        self.setTotalEnthalpy()

    def setTotalEnthalpy(self):

        self.Ht = self.gas.enthalpy_mass() + self.v**2/2.0

    def setMomentumFlux(self):

        self.momFlux = self.gas.pressure() + self.gas.density()*self.v**2

    def setMassFlux(self):

        self.massFlux = self.gas.density()*self.v

    def setVelocity(self, v):

        self.v = v
        self.a = np.sqrt(self.gam*self.R*self.gas.temperature())
        self.M = v/self.a

    def setMach(self, M):

        self.M = M
        self.a = np.sqrt(self.gam*self.R*self.gas.temperature())
        self.v = M*self.a

    def flowSummary(self):
        '''
        Print out a summary of flow properties
        '''

        print "\nrho =", self.gas.density(), "kg/m^3"
        print "T =", self.gas.temperature(), "K"
        print "p =", self.gas.pressure(), "Pa"
        print "R =", self.R, "J/kg/K"
        print "gam =", self.gam

        print "\na =", self.a, "m/s"
        print "v =", self.v, "m/s"
        print "M =", self.M

class shock():

    def __init__(self, inflow=flow(), outflow=flow()):

        self.inflow = inflow
        self.outflow = outflow

        # Set initial guess for p and h
        pinit=self.inflow.gas.density()*self.inflow.v**2.0
        hinit=self.inflow.v**2/2.0
        self.setOutflowState(pinit, hinit)

    def setOutflowState(self, pinit, hinit):
        
        self.outflow.gas.setState_HP(hinit, pinit)
        self.epsilon = self.inflow.gas.density()/self.outflow.gas.density()

    def iterateShock(self):

        p2 = (self.inflow.gas.pressure() + 
              self.inflow.gas.density()*self.inflow.v**2*(1.0-self.epsilon))
        h2 = (self.inflow.gas.enthalpy_mass() +
              0.5*self.inflow.v**2*(1-self.epsilon**2))

        v2 = self.inflow.v*self.epsilon

        self.setOutflowState(p2,h2)
        self.outflow.setVelocity(v2)

    def convergeEpsilon(self,tol=1.0e-5):

        epsilonOld = self.epsilon
        epsilonNew = tol+1.0
        
        while np.abs(epsilonOld-epsilonNew) > tol:

            epsilonOld = self.epsilon
            self.iterateShock()
            epsilonNew = self.epsilon

