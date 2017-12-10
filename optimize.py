#!/usr/bin/python
import sys
sys.path.insert(0, r'./Influenza')
import Optimization

o = Optimization.Optimization()

if not o.options.quiet:
    from miscellany import formatVacSchedule
    print 'Vaccine schedule:\t %s' % formatVacSchedule(o.options.vacSchedule)
    print 'Objective:\t\t %s' % o.options.objective

o.optimize()

if not o.options.quiet:
    o.outputInfo()
