#!/usr/bin/python
import Influenza

s = Influenza.Simulation()
## input parameters 1) vac times, 2) proportion vaccinated across age groups
vacsUsed = s.simulateWithVaccine([0], [[0., 0.05, 0.0364698, 0., 0.]])

print 'Doses:\t\t\t %s' % ', '.join('%g' % v for v in vacsUsed)
print s.outputInfo()
