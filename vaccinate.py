#!/usr/bin/python
#
#
#

import Influenza

s = Influenza.Simulation()

vacsUsed = s.simulateWithVaccine([0], [[0., 0.05, 0.0364698, 0., 0.]], 0.1)

if not s.options.quiet:
    print 'Doses:\t\t\t %s' % ', '.join('%g' % v for v in vacsUsed)
    print

    s.outputInfo()
