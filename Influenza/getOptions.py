#def setRandomSeed(value):
#    # Convert to hex
#    from binascii import hexlify
#    seed = hexlify(value)

#    # Convert into list of ints
#    from sys import maxint
#    S = [long(seed, 16)]
#    while S[-1] > maxint:
#        S.insert(0, int(S[-1] % maxint))
#        S[-1] = S[-1] / maxint
#    S[-1] = int(S[-1])

#    import numpy.random
#    numpy.random.seed(S)

#    return seed

# Set the default RNG seed from the system's random pool
# This is what numpy does by default
from os import urandom
import numpy as np
import random

#defaultSeed = setRandomSeed(urandom(16))
defaultSeed = np.random.seed()

def handleRandomSeed(option, opt_str, value, parser):
    seed = setRandomSeed(value)
    setattr(parser.values, option.dest, seed)

def normalizeObjectiveName(option, opt_str, value, parser):
    from miscellany import normalizeObjectiveName

    setattr(parser.values, option.dest,
            normalizeObjectiveName(value))

def setVacSchedule(option, opt_str, value, parser):
    setattr(parser.values, option.dest, [])
    
    if len(value) == 0:
        parser.print_usage()
        raise ValueError, "Vaccine Information Needed!"

    elif len(value) % 2 != 0:
        parser.print_usage()
        raise ValueError, \
            "Vaccine Information must be in (time, number) pairs!"

    else:
        while len(value) > 0:
            try:
                getattr(parser.values,
                        option.dest).append((float(value.pop(0)),
                                             float(value.pop(0))))
            except ValueError:
                parser.print_usage()
                raise ValueError, \
                    "Vaccine Information must be numbers!"

def getOptions(simType = 'Simulation'):
    import optparse

    p = optparse.OptionParser()

    #p.add_option('-s', '--seed', action = 'callback',
     #            dest = 'seed',
     #            callback = handleRandomSeed,
    #             type = 'string', default = defaultSeed,
    #             help = 'seed for random-number generator')
    p.add_option('-q', '--quiet', action = 'store_true',
                 dest = 'quiet', default = False,
                 help = 'do not output')
    #p.add_option('-d', '--debug', action = 'store_true',
    #             dest = 'debug', default = False,
    #             help = 'output debugging info')

    if simType == 'Optimization':
        p.set_usage(
            'Usage: %prog [options] vacTime0 vacNumber0 [vacTime1 vacNumber1 ...]')
	
        g = optparse.OptionGroup(p, "Optimizitation options")
	
        g.add_option('-o', '--objective', action = 'callback',
                     dest = 'objective',
                     callback = normalizeObjectiveName,
                     type = 'string', default = 'Infections',
                     help = 'objective to optimize over [ default: %default ]')
        g.add_option('-w', '--write', action = 'store_false',
                     dest = 'write', default = True,
                     help =
                     'do not save data from the optimization')
        p.add_option_group(g)

    (options, args) = p.parse_args()
    print ("check!!!"), options,args

    if simType == 'Optimization':
        # vacSchedule is in positional arguments, not options
        # Dummy option container
	print ("vac option"), optparse.Option('-v',dest = 'vacSchedule')
        vacSchedOption = optparse.Option('-v', dest = 'vacSchedule')
	
        setVacSchedule(vacSchedOption, 'vacSchedule', args, p)
    
    return options
