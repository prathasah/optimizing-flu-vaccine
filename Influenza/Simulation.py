import fileIO
import numpy
        
class Simulation:
    def __init__(self, options = None, tMax = 120, paramValues = {}):
        self.tMax = tMax
	

        if options != None:
            self.options = options
	   
        else:
            from getOptions import getOptions
	    

        # Must wait until after options, where RNG seed is set
        Parameters = __import__('Parameters', globals())

        self.parameters = Parameters.Parameters(**paramValues)

        # Initial condition
        self.Y0 = numpy.zeros(8 * self.parameters.ages.size)

        self.hasSolution = False

    def computeR0(self):
        return self.parameters.computeR0()

    def getLastValues(self):
        return (self.SU[-1, :], self.EU[-1, :],
                self.IL[-1, :], self.RU[-1, :],
                self.SV[-1, :], self.EV[-1, :],
                self.IV[-1, :], self.RV[-1, :])

    def updateIC(self):
        if not self.hasSolution:
            # S
	    ## SU
            self.Y0[ 0 : : 8] = \
                     (1 - self.parameters.proportionVaccinated) * self.parameters.population 
	    ## SV
            self.Y0[ 4 : : 8] = \
                     self.parameters.proportionVaccinated * self.parameters.population 

            # E (EU and EV)	    
            self.Y0[ 1 : : 8] = 0.
            self.Y0[ 5 : : 8] = 0.

            # I: Add a single infectious person in each age
            self.Y0[ 2 : : 8] = (1 - self.parameters.proportionVaccinated)
            self.Y0[ 6 : : 8] = self.parameters.proportionVaccinated

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0 : : 8] -= self.Y0[ 2 : : 8]
            self.Y0[ 4 : : 8] -= self.Y0[ 6 : : 8]

            # R (EU and EV)
            self.Y0[ 3 : : 8] = 0.
            self.Y0[ 7 : : 8] = 0.

        else:
            SU, EU, IU, RU, SV, EV, IV, RV = self.getLastValues()
            
            self.Y0[ 0 : : 8] = (1 - self.parameters.proportionVaccinated)* SU
            self.Y0[ 4 : : 8] = SV + self.parameters.proportionVaccinated * SU

            self.Y0[ 1 : : 8] = EU
            self.Y0[ 5 : : 8] = EV
            self.Y0[ 2 : : 8] = IU
            self.Y0[ 6 : : 8] = IV
            self.Y0[ 3 : : 8] = RU
            self.Y0[ 7 : : 8] = RV

    def RHS(self, Y, t):
        '''
        SEIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors
        SU = Y[ 0 : : 8]
        EU = Y[ 1 : : 8]
        IU = Y[ 2 : : 8]
        RU = Y[ 3 : : 8]
        SV = Y[ 4 : : 8]
        EV = Y[ 5 : : 8]
        IV = Y[ 6 : : 8]
        RV = Y[ 7 : : 8]
        
        N = sum(SU + EU + IU + RU + SV + EV + IV + RV)
      
        # The force of infection
        # Lambda[i] = transmissionScaling * susceptibility[i]
        #             * sum(contactMatrix[i, j] * transmissibility[j] * I[j])
        #             / N
        Lambda = self.parameters.transmissionScaling * self.parameters.susceptibility \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility * (IU + IV)) / N
        
        # The right-hand sides
        dSU = - Lambda * SU
        dEU = Lambda * SU - self.parameters.latencyRate * EU
        dIU = self.parameters.latencyRate * EU \
               - (self.parameters.recoveryRate + self.parameters.deathRateU) * IU
        dRU = self.parameters.recoveryRate * IU
        dSV = - (1 - self.vaccineEfficacyVsInfection) * Lambda * SV
        dEV = (1 - self.vaccineEfficacyVsInfection) \
               * Lambda * SV - self.parameters.latencyRate * EV
        dIV = self.parameters.latencyRate * EV \
               - (self.parameters.recoveryRate + self.parameters.deathRateV) * IV
        dRV = self.parameters.recoveryRate * IV
        
        
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 8] = dSU
        dY[ 1 : : 8] = dEU
        dY[ 2 : : 8] = dIU
        dY[ 3 : : 8] = dRU
        dY[ 4 : : 8] = dSV
        dY[ 5 : : 8] = dEV
        dY[ 6 : : 8] = dIV
        dY[ 7 : : 8] = dRV
        
        return dY
    
    def resetSolution(self):
        self.hasSolution = False

    def solve(self, tStart = 0., tEnd = None, tStep = 1.):
        if tEnd == None:
            tEnd = self.tMax

        if self.hasSolution:
            TOld  = self.T.copy()
            SUOld = self.SU.copy()
            EUOld = self.EU.copy()
            IUOld = self.IU.copy()
            RUOld = self.RU.copy()
            SVOld = self.SV.copy()
            EVOld = self.EV.copy()
            IVOld = self.IV.copy()
            RVOld = self.RV.copy()
            
            
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))
        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y = odeint(self.RHS,
                        self.Y0.copy(),
                        self.T,
                        mxstep = 1000)
        Z = self.Y.copy()
        self.SU = Z[:,  0 : : 8]
        self.EU = Z[:,  1 : : 8]
        self.IU = Z[:,  2 : : 8]
        self.RU = Z[:,  3 : : 8]
        self.SV = Z[:,  4 : : 8]
        self.EV = Z[:,  5 : : 8]
        self.IV = Z[:,  6 : : 8]
        self.RV = Z[:,  7 : : 8]

        if self.hasSolution:
            self.T = numpy.hstack((TOld, self.T))

            self.SU = numpy.vstack((SUOld, self.SU))
            self.EU = numpy.vstack((EUOld, self.EU))
            self.IU = numpy.vstack((IUOld, self.IU))
            self.RU = numpy.vstack((RUOld, self.RU))
            self.SV = numpy.vstack((SVOld, self.SV))
            self.EV = numpy.vstack((EVOld, self.EV))
            self.IV = numpy.vstack((IVOld, self.IV))
            self.RV = numpy.vstack((RVOld, self.RV))

        self.hasSolution = True

    def updateStats(self):
        self.NU = self.SU + self.EU + self.IU + self.RU
        self.NV = self.SV + self.EV + self.IV + self.RV
        self.N  = self.NU + self.NV

        self.infectionsU = self.NU[0, :] - self.SU[-1, :]
        self.infectionsV = self.NV[0, :] - self.SV[-1, :]

        # Find duplicate times: these are where vaccination occurs
        for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
            self.infectionsU += self.SU[i + 1, :] - self.SU[i, :]
            self.infectionsV += self.SV[i + 1, :] - self.SV[i, :]

        self.infections  = self.infectionsU + self.infectionsV
        self.totalInfections = self.infections.sum()
        
        self.hospitalizations = self.infections* self.parameters.caseHospitalization
        self.totalHospitalizations = self.hospitalizations.sum()
        
        self.deathsU = self.NU[0, :] - self.NU[-1, :]
        self.deathsV = self.NV[0, :] - self.NV[-1, :]

        # Find duplicate times: these are where vaccination occurs
        for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
            self.deathsU += self.NU[i + 1, :] - self.NU[i, :]
            self.deathsV += self.NV[i + 1, :] - self.NV[i, :]

        self.deaths  = self.deathsU + self.deathsV
        self.totalDeaths = self.deaths.sum()
        
        
    def simulate(self):
        self.updateIC()
        self.solve()
        self.updateStats()

    def updateProportionVaccinated(self, PVPWVal):
        # Update propotion vaccinated
        self.parameters.proportionVaccinatedPW.values = PVPWVal
	## extend to full ages groups. Proportions calculated by multiplying PVPWVal 
	##values with the matrix defined in S.130
        self.parameters.proportionVaccinated = self.parameters.proportionVaccinatedPW.full(self.parameters.ages)
	if self.hasSolution:
            IC = self.getLastValues()
            vacsUsed = (self.parameters.proportionVaccinated * IC[0]).sum()

        else:
            vacsUsed = (self.parameters.proportionVaccinated
                        * self.parameters.population).sum() 
	

        # Update initial condition for ODEs
        self.updateIC()

        return vacsUsed
        
    def simulateWithVaccine(self, vacTimes, PVPWVals, vacEfficacy):
	
	self.vaccineEfficacyVsInfection = vacEfficacy * self.parameters.relative_vaccineEfficacyVsInfection
	
        nVacRounds = len(vacTimes)

        # Convert flat vector to 2-D array
        if numpy.ndim(PVPWVals) != 2:
            PVPWVals = numpy.asarray(PVPWVals).reshape(
                (nVacRounds,
                 self.parameters.proportionVaccinatedLength))

        self.resetSolution()

        if not (vacTimes[0] == 0.):
            # Solve from 0 to time of first vaccine
            # Set to no vaccine
            self.updateProportionVaccinated(
                [0.] * self.parameters.proportionVaccinatedLength)
            self.solve(tEnd = vacTimes[0])

        vacsUsed = numpy.empty(nVacRounds)
	for vacRound in range(nVacRounds):
            # Vaccinate the population
            vacsUsed[vacRound] = self.updateProportionVaccinated(
                PVPWVals[vacRound])

            # Run from now until
            if vacRound + 1 < nVacRounds:
                # Next vaccine
                tEnd = vacTimes[vacRound + 1]
            else:
                # End time
                tEnd = self.tMax
	    self.solve(tStart = vacTimes[vacRound], tEnd = tEnd)

        self.updateStats()

        return vacsUsed

    def outputSolution(self):
        for (i, t) in enumerate(self.T):
            print '%g' % t,
            print '\t'.join(map(lambda f: ('%g' % f),
                                (sum(self.EU[i, :]),
                                 sum(self.IU[i, :]),
                                 sum(self.RU[i, :]),
                                 sum(self.SV[i, :]),
                                 sum(self.EV[i, :]),
                                 sum(self.IV[i, :]),
                                 sum(self.RV[i, :]))))
    
    def outputInfo(self):
        print 'R0:\t\t\t %g' % self.parameters.R0
        print ('Infections (in millions):'),  self.totalInfections/1000000.
        print 'Deaths:\t\t\t %g' % self.totalDeaths
        print 'Hospitalizations:\t %g' % self.totalHospitalizations
        #print ('Age-specific infections:'), list(self.infections)
	print ("unvaccinated"),self.infectionsU.sum()
	print ("vaccinated"),self.infectionsV.sum()
	print ("total pop size"), self.parameters.population.sum()
	print ("total unvaccinated"),  ((1 - self.parameters.proportionVaccinated) * self.parameters.population).sum()
	print ("total vaccinated"),  ((self.parameters.proportionVaccinated) * self.parameters.population).sum()


    def short_output(self):
	return list(self.infections), list(self.hospitalizations)

    def debug_info(self):
	#print ("recovery rate ="), self.parameters.recoveryRate 
	#print ("latency rate=="), self.parameters.latencyRate
	return self.parameters.deathRateU
	

    def vaccinated_output(self):
        return self.parameters.population.sum(),((1 - self.parameters.proportionVaccinated) * self.parameters.population).sum(), ((self.parameters.proportionVaccinated) * self.parameters.population).sum(), self.infectionsU.sum(), self.infectionsV.sum()

    def getDumpData(self, runData = True):
        dumpData = fileIO.dumpContainer()

        dumpData.settings = fileIO.dumpContainer()
        dumpData.settings.seed = self.options.seed
        dumpData.settings.tMax = self.tMax

        dumpData.parameters = self.parameters.getDumpData()

        dumpData.stats = fileIO.dumpContainer()
        dumpData.stats.infectionsUL = self.infectionsUL
        dumpData.stats.infectionsVL = self.infectionsVL
        dumpData.stats.infectionsUH = self.infectionsUH
        dumpData.stats.infectionsVH = self.infectionsVH
        dumpData.stats.infectionsU = self.infectionsU
        dumpData.stats.infectionsV = self.infectionsV
        dumpData.stats.infectionsL = self.infectionsL
        dumpData.stats.infectionsH = self.infectionsH
        dumpData.stats.infections = self.infections
        dumpData.stats.totalInfections = self.totalInfections
        dumpData.stats.deathsUL = self.deathsUL
        dumpData.stats.deathsVL = self.deathsVL
        dumpData.stats.deathsUH = self.deathsUH
        dumpData.stats.deathsVH = self.deathsVH
        dumpData.stats.deathsU = self.deathsU
        dumpData.stats.deathsV = self.deathsV
        dumpData.stats.deathsL = self.deathsL
        dumpData.stats.deathsH = self.deathsH
        dumpData.stats.deaths = self.deaths
        dumpData.stats.totalDeaths = self.totalDeaths
        dumpData.stats.hospitalizationsL = self.hospitalizationsL
        dumpData.stats.hospitalizationsH = self.hospitalizationsH
        dumpData.stats.hospitalizations = self.hospitalizations
        dumpData.stats.totalHospitalizations = self.totalHospitalizations
        dumpData.stats.yearsOfLifeLost = self.YLL
        dumpData.stats.contingentValuation = self.CV
        #dumpData.stats.dollarCost = self.totalCost

        if runData:
            dumpData.runData = fileIO.dumpContainer()
            for key in ('T',
                        'SUL', 'SVL', 'SUH', 'SVH',
                        'EUL', 'EVL', 'EUH', 'EVH',
                        'IUL', 'IVL', 'IUH', 'IVH',
                        'RUL', 'RVL', 'RUH', 'RVH',
                        'NUL', 'NVL', 'NUH', 'NVH',
                        'NU', 'NV', 'NL', 'NH', 'N'):
                setattr(dumpData.runData, key, getattr(self, key))

        return dumpData

    dump = fileIO.dump
    dumps = fileIO.dumps
    
    def plot(self, show = True):
        import pylab
        
        pylab.plot(self.T,
                   (self.IUL + self.IVL + self.IUH + self.IVH).sum(1)
                   / self.N.sum(1))
        pylab.xlabel('Time (days)')
        pylab.ylabel('Proportion infectious')

        pylab.ylim(ymin = 0)

        if show:
            pylab.show()
