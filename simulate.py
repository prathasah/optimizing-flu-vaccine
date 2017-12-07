#!/usr/bin/python
#
#
#

import Influenza

s = Influenza.Simulation()

s.simulate()

if not s.options.quiet:
    s.outputInfo()
