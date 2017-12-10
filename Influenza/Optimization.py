import fileIO

import numpy

class Optimization:
    objectiveMap = {'Infections': 'totalInfections',
                    'Deaths': 'totalDeaths',
                    'Hospitalizations': 'totalHospitalizations',
                    'YLL': 'YLL',
                    'Contingent': 'CV',
                    'Cost': 'totalCost'}

    def __init__(self, options = None, optimRuns = 3,
                 *args, **kwargs):
        self.optimRuns = optimRuns

        if options != None:
            self.options = options
        else:
            from getOptions import getOptions
            self.options = getOptions('Optimization')
        
        self.vacTimes = numpy.array([v[0] for v in self.options.vacSchedule])
        self.vacNumbers = numpy.array([v[1] for v in self.options.vacSchedule])
	self.vacEfficacy = numpy.array([v[2] for v in self.options.vacSchedule])
        self.nVacRounds = len(self.vacTimes)

        self.objective = self.options.objective

        from Simulation import run_Simulation

        self.s = run_Simulation(options = self.options, *args, **kwargs)

        self.PVUsed = None

        self.vacsUsed = numpy.array([0] * self.nVacRounds)

    def solve(self, PVPWVals):
        # Only update for new PVPWVals
        if numpy.any(PVPWVals != self.PVUsed):
	    self.vacsUsed = self.s.simulateWithVaccine(self.vacTimes,
                                                       PVPWVals, self.vacEfficacy[0])

            self.PVUsed = PVPWVals.copy()

    def evaluateObjective(self, PVPWVals):
        self.solve(PVPWVals)
        return getattr(self.s, self.objectiveMap[self.objective])

    def totalVacsConditions(self, PVPWVals):
        self.solve(PVPWVals)
        return (self.vacNumbers - self.vacsUsed)

    def totalVacsCondition(self, i):
        return lambda PVPWVals: self.totalVacsConditions(PVPWVals)[i]

    def lowerCondition(self, i):
        return lambda PVPWVals: PVPWVals[i]

    def upperCondition(self, i):
        return lambda PVPWVals: 1.0 - PVPWVals[i]

    def optimize(self):
        from scipy.optimize import fmin_cobyla

        conds = [self.totalVacsCondition(i) for i in
                 range(self.nVacRounds)]
        conds.extend([self.lowerCondition(i) for i in
                      range(self.s.parameters.proportionVaccinatedLength
                            * self.nVacRounds)])
        conds.extend([self.upperCondition(i) for i in
                      range(self.s.parameters.proportionVaccinatedLength
                            * self.nVacRounds)])
        
        minObjective = None

        for i in range(self.optimRuns):
            if self.options.debug:
                print 'Run %d' % (i + 1)

            PV0 = 0.5 * numpy.random.rand(
                self.s.parameters.proportionVaccinatedLength
                * self.nVacRounds)

            PVPWValsOpt = fmin_cobyla(self.evaluateObjective,
                                      PV0,
                                      conds,
                                      maxfun = 10000,
                                      rhobeg = 0.5,
                                      iprint = 0)

            if (minObjective == None) \
                    or (self.evaluateObjective(PVPWValsOpt) < minObjective):
                if self.options.debug:
                    print 'Optimum improved.'
                
                minObjective = self.evaluateObjective(PVPWValsOpt)
                self.PVBest = PVPWValsOpt
        
	self.solve(self.PVBest)

        if self.options.write:
            self.dump(self.openDumpFile())

    def outputInfo(self):
        print '\nValues:\t\t',
        for PVI in self.PVBest:
            print ' %g' % PVI,
        print '\n'
        
        print 'Doses:\t\t\t %s' % ', '.join(['%g' % v for v in self.vacsUsed])
        print 'Objective:\t\t %g' % self.evaluateObjective(self.PVBest)
        
        print
        
        self.s.outputInfo()

    def getDumpData(self):
        dumpData = fileIO.dumpContainer()

        dumpData.settings = fileIO.dumpContainer()
        dumpData.settings.objective = self.objective
        dumpData.settings.seed = self.options.seed
        dumpData.settings.vacTimes = self.vacTimes
        dumpData.settings.vacNumbers = self.vacNumbers
        dumpData.settings.nVacRounds = self.nVacRounds
        dumpData.settings.optimRuns = self.optimRuns

        # Optimal vaccination
        dumpData.proportionVaccinated = self.PVBest
        dumpData.objectiveValue = self.evaluateObjective(self.PVBest)
        dumpData.vacsUsed = self.vacsUsed

        dumpData.simulation = self.s.getDumpData(runData = False)

        return dumpData

    dump = fileIO.dump
    dumps = fileIO.dumps
    
    def encodeVacSchedule(self):
        return '_'.join(['%g,%d' % (t, int(v)) for (t, v)
                         in self.options.vacSchedule])

    def getDumpDir(self):
        from os.path import join

        return join('data', '%s_%s'
                    % (self.encodeVacSchedule(), self.objective))


    def openDumpFile(self):
        import os.path

        dumpDir = self.getDumpDir()

        if not os.path.exists(dumpDir):
            from os import makedirs
            makedirs(dumpDir, 0755)

        fileName = os.path.join(dumpDir, '%s' % self.options.seed)

        if self.options.debug:
            print 'Writing output to file %s' % fileName

        if os.path.exists(fileName):
            print 'Overwriting output file %s!' % fileName

        return open(fileName, 'wb')
