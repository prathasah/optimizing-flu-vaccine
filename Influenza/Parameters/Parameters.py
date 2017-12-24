from PiecewiseAgeParameter import PiecewiseAgeParameter, PiecewiseAgeRate
from ages import ages, vaccinationAges
import demography
import epidemiology
import costs
import os
import numpy as np
#from .. import fileIO


import numpy

import UserDict

class ParamDict(UserDict.UserDict):
    def valueOrAttrFromOther(self, key, other):
        '''
        Return value if key is in paramValues dict or return
        default as attribute from object other.
        '''
	
	reload(epidemiology)
	reload(costs)
	return getattr(other, key)

class Parameters:
    def setPWAttr(self, namePW, value):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name.
        '''
        assert isinstance(value, PiecewiseAgeParameter)
        assert namePW.endswith('PW')
        name = namePW[ : -2]
        setattr(self, namePW, value)
        setattr(self, name, value.full(self.ages))
    
    def setPWAttrFromPassedOrOther(self, other, namePW):
        '''
        Set a piecewise attribute as self.namePW
        and its expanded value as self.name,
        taking values from passed paramValues
        or from attributes of object other.
        '''
	self.setPWAttr(namePW,self.passedParamValues.valueOrAttrFromOther(namePW, other))
	

    def setAttrFromPassedOrOther(self, other, name):
        '''
        Set an attribute as self.name,
        taking value from passed paramValues
        or from attributes of object other.
        '''
	setattr(self, name, 
                self.passedParamValues.valueOrAttrFromOther(name, other))

##################################################################

    def __init__(self, **paramValues):

	self.passedParamValues = ParamDict(paramValues)	

	self.ages = numpy.array(ages)

        # Load in parameters and expand as necessary
	# Go thrrough each files
        for m in (demography, epidemiology, costs):
	    # list all the modules
            for p in dir(m):
		#if module returns a numbers, then..
                if isinstance(getattr(m, p),(float, int)):
		    self.setAttrFromPassedOrOther(m, p)
		##if it is an agespecific parameter then..
                elif isinstance(getattr(m, p),
                                PiecewiseAgeParameter): 
		    self.setPWAttrFromPassedOrOther(m, p)

	self.vaccineEfficacyVsInfection = self.passedParamValues["vacEfficacy"] * self.relative_vaccineEfficacyVsInfection

	self.vaccineCoverage = [(self.passedParamValues["vacCoverage"]*num/np.mean(self.relative_coverage)) for num in self.relative_coverage]


        # Compute mortality rates
        # from case mortalities

	# Mortality of *vaccinated* risk individuals
        self.caseMortalityV = self.caseMortality \
                               * (1 - self.vaccineEfficacyVsDeath)
	
	##Death rate of unvaccinated individuals
        self.deathRateU = (self.recoveryRate * self.caseMortality) / (1 - self.caseMortality)
	##Death rate of vaccinated individuals
        self.deathRateV = (self.recoveryRate * self.caseMortalityV)/ (1 - self.caseMortalityV)
	
        # Set up proportion vaccinated vectors
        self.proportionVaccinatedPW = PiecewiseAgeRate(
            [0.0] * len(vaccinationAges),
            vaccinationAges)
	self.proportionVaccinated = self.proportionVaccinatedPW.full(self.ages)
	self.proportionVaccinatedLength = len(vaccinationAges)
        

        # Get contact matrix
        if 'contactMatrix' in paramValues:
            self.contactMatrix = paramValues.get('contactMatrix')
        else:
            from sys import modules
            import os.path
            import cPickle
            modulePath = os.path.dirname(modules[self.__module__].__file__)
            contactMatrixFile = os.path.join(modulePath, 'contactMatrix.p')
            self.contactMatrix = cPickle.load(open(contactMatrixFile))

        # One last parameter to fit to R0
        self.transmissionScaling = 1.0
	self.transmissionScaling *=  self.R0 / self.computeR0()
	
        
    def computeR0(self):
	#normalized population size for each age groups
        s0 = self.population / sum(self.population)
	sU0 = s0 * (1 - self.proportionVaccinated)
	sV0 = s0 * self.proportionVaccinated
	

        FU = self.transmissionScaling \
              * numpy.outer(self.susceptibility * sU0,
                            self.transmissibility) * self.contactMatrix
	FV = self.transmissionScaling \
              * numpy.outer((1 - self.vaccineEfficacyVsInfection)
                            * self.susceptibility * sV0,
                            self.transmissibility) * self.contactMatrix
        
        F = numpy.vstack((numpy.hstack((FU, FU)),
                          numpy.hstack((FV, FV))))

        V = numpy.diag(numpy.hstack(
            (self.recoveryRate + self.deathRateU,
             self.recoveryRate + self.deathRateV)))

        G = numpy.dot(F, numpy.linalg.inv(V))

        (Lambda, Nu) = numpy.linalg.eig(G)

        i = numpy.argmax(Lambda.real)
        LambdaMax = numpy.real_if_close(Lambda[i])
        
        assert numpy.isreal(LambdaMax), \
               "Complex maximal eigenvalue %s!" % LambdaMax

        return LambdaMax

    def getDumpData(self, runData = True):
        dumpData = fileIO.dumpContainer()

        dumpData.ages = self.ages
        dumpData.population = self.population
        dumpData.R0 = self.R0
        dumpData.recoveryRate = self.recoveryRatePW
        dumpData.latencyRate = self.latencyRatePW
        dumpData.susceptibility = self.susceptibilityPW
        dumpData.transmissibility = self.transmissibilityPW
        dumpData.relative_vaccineEfficacyVsInfection = \
            self.relative_vaccineEfficacyVsInfectionPW
        dumpData.vaccineEfficacyVsDeath = \
            self.vaccineEfficacyVsDeathPW
        dumpData.caseMortality = self.caseMortalityPW
        dumpData.highRiskRelativeCaseMortality = \
            self.highRiskRelativeCaseMortalityPW
        dumpData.contactMatrix = self.contactMatrix
        dumpData.expectationOfLife = self.expectationOfLife
        dumpData.contingentValue = self.contingentValue
        dumpData.caseHospitalization = \
            self.caseHospitalizationPW
        dumpData.highRiskRelativeCaseHospitalization = \
            self.highRiskRelativeCaseHospitalizationPW

        return dumpData

    #dump = fileIO.dump
