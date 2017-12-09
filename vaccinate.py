#!/usr/bin/python
#
#
#

import Influenza

s = Influenza.Simulation()

vacsUsed = s.simulateWithVaccine([0], [[0.7012, 0.5514, 0.3264, 0.4528, 0.6532]], 0.1)

#if not s.options.quiet:
#    print 'Doses:\t\t\t %s' % ', '.join('%g' % v for v in vacsUsed)
#    print

s.outputInfo()
