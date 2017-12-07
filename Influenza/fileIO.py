import cPickle

class dumpContainer:
    """Empty container"""
    pass

def dump(obj, fileDescriptor):
    cPickle.dump(obj.getDumpData(), fileDescriptor)

def dumps(obj):
    return cPickle.dumps(obj.getDumpData())

load = cPickle.load

loads = cPickle.loads
