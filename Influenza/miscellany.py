def makeXTicks(ages):
    ticks = []
    
    for i in range(len(ages) - 1):
        begin = ages[i]
        end = ages[i + 1] - 1
        
        if end > begin:
            ticks.append(str(begin) + '-' + str(end))
        else:
            ticks.append(str(begin))
    
    ticks.append(str(ages[-1]) + '+')
    
    return ticks

def normalizeObjectiveName(objective):
    if objective.strip().lower() == 'yll':
        return objective.strip().upper()
    else:
        return objective.strip().capitalize()

def formatVacSchedule(vacSchedule):
    s = []
    
    for (t, v, e) in vacSchedule:
        s.append('(%g, %g, %g)' % (t, v, e))

    return ', '.join(s)
