import fileIO

import numpy
        
class Simulation:
    def __init__(self, options = None, tMax = 1000, paramValues = {}):
        self.tMax = tMax

        if options != None:
            self.options = options
        else:
            from getOptions import getOptions
            self.options = getOptions('Simulation')

        # Must wait until after options, where RNG seed is set
        Parameters = __import__('Parameters', globals())

        self.parameters = Parameters.Parameters(**paramValues)

        # Initial condition
        self.Y0 = numpy.zeros(16 * self.parameters.ages.size)

        self.hasSolution = False

    def computeR0(self):
        return self.parameters.computeR0()

    def getLastValues(self):
        return (self.SUL[-1, :], self.EUL[-1, :],
                self.IUL[-1, :], self.RUL[-1, :],
                self.SVL[-1, :], self.EVL[-1, :],
                self.IVL[-1, :], self.RVL[-1, :],
                self.SUH[-1, :], self.EUH[-1, :],
                self.IUH[-1, :], self.RUH[-1, :],
                self.SVH[-1, :], self.EVH[-1, :],
                self.IVH[-1, :], self.RVH[-1, :])

    def updateIC(self):
        if not self.hasSolution:
            # S
            self.Y0[ 0 : : 16] = \
                     (1 - self.parameters.proportionVaccinatedL) \
                     * self.parameters.population \
                     * (1 - self.parameters.proportionHighRisk)
            self.Y0[ 4 : : 16] = \
                     self.parameters.proportionVaccinatedL \
                     * self.parameters.population \
                     * (1 - self.parameters.proportionHighRisk)
            self.Y0[ 8 : : 16] = \
                     (1 - self.parameters.proportionVaccinatedH) \
                     * self.parameters.population \
                     * self.parameters.proportionHighRisk
            self.Y0[12 : : 16] = \
                       self.parameters.proportionVaccinatedH \
                       * self.parameters.population \
                       * self.parameters.proportionHighRisk

            # E
            self.Y0[ 1 : : 16] = 0.
            self.Y0[ 5 : : 16] = 0.
            self.Y0[ 9 : : 16] = 0.
            self.Y0[13 : : 16] = 0.

            # I: Add a single infectious person in each age
            self.Y0[ 2 : : 16] = \
                     (1 - self.parameters.proportionVaccinatedL)
            self.Y0[ 6 : : 16] = self.parameters.proportionVaccinatedL
            self.Y0[10 : : 16] = \
                       (1 - self.parameters.proportionVaccinatedH)
            self.Y0[14 : : 16] = self.parameters.proportionVaccinatedH

            # S: Remove those new infectious people from the susceptibles
            self.Y0[ 0 : : 16] -= self.Y0[ 2 : : 16]
            self.Y0[ 4 : : 16] -= self.Y0[ 6 : : 16]
            self.Y0[ 8 : : 16] -= self.Y0[10 : : 16]
            self.Y0[12 : : 16] -= self.Y0[14 : : 16]

            # R
            self.Y0[ 3 : : 16] = 0.
            self.Y0[ 7 : : 16] = 0.
            self.Y0[11 : : 16] = 0.
            self.Y0[15 : : 16] = 0.

        else:
            SUL, EUL, IUL, RUL, SVL, EVL, IVL, RVL, \
                SUH, EUH, IUH, RUH, SVH, EVH, IVH, RVH = \
                self.getLastValues()
            
            self.Y0[ 0 : : 16] = (1 - self.parameters.proportionVaccinatedL) \
                                 * SUL
            self.Y0[ 4 : : 16] = SVL \
                                 + self.parameters.proportionVaccinatedL * SUL
            self.Y0[ 8 : : 16] = (1 - self.parameters.proportionVaccinatedH) \
                                 * SUH
            self.Y0[12 : : 16] = SVH \
                                 + self.parameters.proportionVaccinatedH * SUH

            self.Y0[ 1 : : 16] = EUL
            self.Y0[ 5 : : 16] = EVL
            self.Y0[ 9 : : 16] = EUH
            self.Y0[13 : : 16] = EVH
            self.Y0[ 2 : : 16] = IUL
            self.Y0[ 6 : : 16] = IVL
            self.Y0[10 : : 16] = IUH
            self.Y0[14 : : 16] = IVH
            self.Y0[ 3 : : 16] = RUL
            self.Y0[ 7 : : 16] = RVL
            self.Y0[11 : : 16] = RUH
            self.Y0[15 : : 16] = RVH

    def RHS(self, Y, t):
        '''
        SEIR model with multiple host types.
        
        This function gives the right-hand sides of the ODEs.
        '''
        
        # Convert vector to meaningful component vectors
        SUL = Y[ 0 : : 16]
        EUL = Y[ 1 : : 16]
        IUL = Y[ 2 : : 16]
        RUL = Y[ 3 : : 16]
        SVL = Y[ 4 : : 16]
        EVL = Y[ 5 : : 16]
        IVL = Y[ 6 : : 16]
        RVL = Y[ 7 : : 16]
        SUH = Y[ 8 : : 16]
        EUH = Y[ 9 : : 16]
        IUH = Y[10 : : 16]
        RUH = Y[11 : : 16]
        SVH = Y[12 : : 16]
        EVH = Y[13 : : 16]
        IVH = Y[14 : : 16]
        RVH = Y[15 : : 16]
        N = sum(SUL + EUL + IUL + RUL + SVL + EVL + IVL + RVL
                + SUH + EUH + IUH + RUH + SVH + EVH + IVH + RVH)
      
        # The force of infection
        # Lambda[i] = transmissionScaling * susceptibility[i]
        #             * sum(contactMatrix[i, j] * transmissibility[j] * I[j])
        #             / N
        Lambda = self.parameters.transmissionScaling \
                 * self.parameters.susceptibility \
                 * numpy.dot(self.parameters.contactMatrix,
                             self.parameters.transmissibility \
                             * (IUL + IVL + IUH + IVH)) / N
        
        # The right-hand sides
        dSUL = - Lambda * SUL
        dEUL = Lambda * SUL - self.parameters.latencyRate * EUL
        dIUL = self.parameters.latencyRate * EUL \
               - (self.parameters.recoveryRate + self.parameters.deathRateUL) \
               * IUL
        dRUL = self.parameters.recoveryRate * IUL
        dSVL = - (1 - self.parameters.relative_vaccineEfficacyVsInfection) \
               * Lambda * SVL
        dEVL = (1 - self.parameters.relative_vaccineEfficacyVsInfection) \
               * Lambda * SVL \
               - self.parameters.latencyRate * EVL
        dIVL = self.parameters.latencyRate * EVL \
               - (self.parameters.recoveryRate + self.parameters.deathRateVL) \
               * IVL
        dRVL = self.parameters.recoveryRate * IVL
        dSUH = - Lambda * SUH
        dEUH = Lambda * SUH - self.parameters.latencyRate * EUH
        dIUH = self.parameters.latencyRate * EUH \
               - (self.parameters.recoveryRate + self.parameters.deathRateUH) \
               * IUH
        dRUH = self.parameters.recoveryRate * IUH
        dSVH = - (1 - self.parameters.relative_vaccineEfficacyVsInfection) \
               * Lambda * SVH
        dEVH = (1 - self.parameters.relative_vaccineEfficacyVsInfection) \
               * Lambda * SVH \
               - self.parameters.latencyRate * EVH
        dIVH = self.parameters.latencyRate * EVH \
               - (self.parameters.recoveryRate + self.parameters.deathRateVH) \
               * IVH
        dRVH = self.parameters.recoveryRate * IVH
        
        # Convert meaningful component vectors into a single vector
        dY = numpy.empty(Y.size, dtype = float)
        dY[ 0 : : 16] = dSUL
        dY[ 1 : : 16] = dEUL
        dY[ 2 : : 16] = dIUL
        dY[ 3 : : 16] = dRUL
        dY[ 4 : : 16] = dSVL
        dY[ 5 : : 16] = dEVL
        dY[ 6 : : 16] = dIVL
        dY[ 7 : : 16] = dRVL
        dY[ 8 : : 16] = dSUH
        dY[ 9 : : 16] = dEUH
        dY[10 : : 16] = dIUH
        dY[11 : : 16] = dRUH
        dY[12 : : 16] = dSVH
        dY[13 : : 16] = dEVH
        dY[14 : : 16] = dIVH
        dY[15 : : 16] = dRVH
        
        return dY
    
    def resetSolution(self):
        self.hasSolution = False

    def solve(self, tStart = 0., tEnd = None, tStep = 1.):
        if tEnd == None:
            tEnd = self.tMax

        if self.hasSolution:
            TOld  = self.T.copy()
            SULOld = self.SUL.copy()
            EULOld = self.EUL.copy()
            IULOld = self.IUL.copy()
            RULOld = self.RUL.copy()
            SVLOld = self.SVL.copy()
            EVLOld = self.EVL.copy()
            IVLOld = self.IVL.copy()
            RVLOld = self.RVL.copy()
            SUHOld = self.SUH.copy()
            EUHOld = self.EUH.copy()
            IUHOld = self.IUH.copy()
            RUHOld = self.RUH.copy()
            SVHOld = self.SVH.copy()
            EVHOld = self.EVH.copy()
            IVHOld = self.IVH.copy()
            RVHOld = self.RVH.copy()
            
        # Time vector for solution
        self.T = numpy.hstack((numpy.arange(tStart, tEnd, tStep), tEnd))
        
        # Integrate the ODE
        from scipy.integrate import odeint
        self.Y = odeint(self.RHS,
                        self.Y0.copy(),
                        self.T,
                        mxstep = 1000)
        
        Z = self.Y.copy()
        self.SUL = Z[:,  0 : : 16]
        self.EUL = Z[:,  1 : : 16]
        self.IUL = Z[:,  2 : : 16]
        self.RUL = Z[:,  3 : : 16]
        self.SVL = Z[:,  4 : : 16]
        self.EVL = Z[:,  5 : : 16]
        self.IVL = Z[:,  6 : : 16]
        self.RVL = Z[:,  7 : : 16]
        self.SUH = Z[:,  8 : : 16]
        self.EUH = Z[:,  9 : : 16]
        self.IUH = Z[:, 10 : : 16]
        self.RUH = Z[:, 11 : : 16]
        self.SVH = Z[:, 12 : : 16]
        self.EVH = Z[:, 13 : : 16]
        self.IVH = Z[:, 14 : : 16]
        self.RVH = Z[:, 15 : : 16]

        if self.hasSolution:
            self.T = numpy.hstack((TOld, self.T))

            self.SUL = numpy.vstack((SULOld, self.SUL))
            self.EUL = numpy.vstack((EULOld, self.EUL))
            self.IUL = numpy.vstack((IULOld, self.IUL))
            self.RUL = numpy.vstack((RULOld, self.RUL))
            self.SVL = numpy.vstack((SVLOld, self.SVL))
            self.EVL = numpy.vstack((EVLOld, self.EVL))
            self.IVL = numpy.vstack((IVLOld, self.IVL))
            self.RVL = numpy.vstack((RVLOld, self.RVL))
            self.SUH = numpy.vstack((SUHOld, self.SUH))
            self.EUH = numpy.vstack((EUHOld, self.EUH))
            self.IUH = numpy.vstack((IUHOld, self.IUH))
            self.RUH = numpy.vstack((RUHOld, self.RUH))
            self.SVH = numpy.vstack((SVHOld, self.SVH))
            self.EVH = numpy.vstack((EVHOld, self.EVH))
            self.IVH = numpy.vstack((IVHOld, self.IVH))
            self.RVH = numpy.vstack((RVHOld, self.RVH))

        self.hasSolution = True

    def updateStats(self):
        self.NUL = self.SUL + self.EUL + self.IUL + self.RUL
        self.NVL = self.SVL + self.EVL + self.IVL + self.RVL
        self.NL  = self.NUL + self.NVL
        self.NUH = self.SUH + self.EUH + self.IUH + self.RUH
        self.NVH = self.SVH + self.EVH + self.IVH + self.RVH
        self.NH  = self.NUH + self.NVH
        self.NU  = self.NUL + self.NUH
        self.NV  = self.NVL + self.NVH
        self.N   = self.NL  + self.NH

        self.infectionsUL = self.NUL[0, :] - self.SUL[-1, :]
        self.infectionsVL = self.NVL[0, :] - self.SVL[-1, :]
        self.infectionsUH = self.NUH[0, :] - self.SUH[-1, :]
        self.infectionsVH = self.NVH[0, :] - self.SVH[-1, :]

        # Find duplicate times: these are where vaccination occurs
        for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
            self.infectionsUL += self.SUL[i + 1, :] - self.SUL[i, :]
            self.infectionsVL += self.SVL[i + 1, :] - self.SVL[i, :]
            self.infectionsUH += self.SUH[i + 1, :] - self.SUH[i, :]
            self.infectionsVH += self.SVH[i + 1, :] - self.SVH[i, :]

        self.infectionsL  = self.infectionsUL + self.infectionsVL
        self.infectionsH  = self.infectionsUH + self.infectionsVH
        self.infectionsU  = self.infectionsUL + self.infectionsUH
        self.infectionsV  = self.infectionsVL + self.infectionsVH
        self.infections   = self.infectionsL  + self.infectionsH
	self.totalInfections = self.infections.sum()
        
        self.hospitalizationsL = self.infectionsL \
                                 * self.parameters.caseHospitalizationL
        self.hospitalizationsH = self.infectionsH \
                                 * self.parameters.caseHospitalizationH
        self.hospitalizations  = self.hospitalizationsL \
                                 + self.hospitalizationsH
        self.totalHospitalizations = self.hospitalizations.sum()
        
        self.deathsUL = self.NUL[0, :] - self.NUL[-1, :]
        self.deathsVL = self.NVL[0, :] - self.NVL[-1, :]
        self.deathsUH = self.NUH[0, :] - self.NUH[-1, :]
        self.deathsVH = self.NVH[0, :] - self.NVH[-1, :]

        # Find duplicate times: these are where vaccination occurs
        for i in numpy.arange(len(self.T)).compress(numpy.diff(self.T) == 0):
            # Update for vaccinated people
            self.deathsUL += self.NUL[i + 1, :] - self.NUL[i, :]
            self.deathsVL += self.NVL[i + 1, :] - self.NVL[i, :]
            self.deathsUH += self.NUH[i + 1, :] - self.NUH[i, :]
            self.deathsVH += self.NVH[i + 1, :] - self.NVH[i, :]

        self.deathsL  = self.deathsUL + self.deathsVL
        self.deathsH  = self.deathsUH + self.deathsVH
        self.deathsU  = self.deathsUL + self.deathsUH
        self.deathsV  = self.deathsVL + self.deathsVH
        self.deaths   = self.deathsL  + self.deathsH
        self.totalDeaths = self.deaths.sum()
        
        self.YLL = numpy.dot(self.parameters.expectationOfLife,
                             self.deaths)
        
        self.CV = numpy.dot(self.parameters.contingentValue,
                            self.deaths)
        
    def simulate(self):
        self.updateIC()
        self.solve()
        self.updateStats()

    def updateProportionVaccinated(self, PVPWVal):
        # Update propotion vaccinated
        self.parameters.proportionVaccinatedLPW.values = \
           PVPWVal[: self.parameters.proportionVaccinatedLLength]
	## extend to full ages groups. Proportions calculated by multiplying PVPWVal first five
	##values with the matrix defined in S.130
        self.parameters.proportionVaccinatedL = \
           self.parameters.proportionVaccinatedLPW.full(self.parameters.ages)
	
	if self.hasSolution:
            IC = self.getLastValues()
            vacsUsed = (self.parameters.proportionVaccinatedL * IC[0]).sum()

        else:
            vacsUsed = (self.parameters.proportionVaccinatedL
                        * self.parameters.population).sum() 

        # Update initial condition for ODEs
        self.updateIC()

        return vacsUsed
        
    def simulateWithVaccine(self, vacTimes, PVPWVals):
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
                                (sum(self.EUL[i, :]),
                                 sum(self.IUL[i, :]),
                                 sum(self.RUL[i, :]),
                                 sum(self.SVL[i, :]),
                                 sum(self.EVL[i, :]),
                                 sum(self.IVL[i, :]),
                                 sum(self.RVL[i, :]),
                                 sum(self.SUH[i, :]),
                                 sum(self.EUH[i, :]),
                                 sum(self.IUH[i, :]),
                                 sum(self.RUH[i, :]),
                                 sum(self.SVH[i, :]),
                                 sum(self.EVH[i, :]),
                                 sum(self.IVH[i, :]),
                                 sum(self.RVH[i, :]))))
    
    def outputInfo(self):
        print 'R0:\t\t\t %g' % self.parameters.R0
        print 'Infections:\t\t %g' % self.totalInfections
        print 'Deaths:\t\t\t %g' % self.totalDeaths
        print 'Hospitalizations:\t %g' % self.totalHospitalizations
        print 'YLL:\t\t\t %g' % self.YLL
        print 'Contingent:\t\t %g' % self.CV
	print ('Age-specific infections:'), list(self.infections), sum(list(self.infections))

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
