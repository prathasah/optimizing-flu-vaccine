import numpy


class PiecewiseAgeParameter:
    """
    Has values[0] between ages[0] and ages[1],
        values[1] between ages[1] and ages[2],
        ...
    """
    def __init__(self, values, ages):
        assert (len(values) == len(ages)), \
               "Values and ages are different lengths!"

        self.values = values
        self.ages = ages
    
    def makeTransMat(self, fullAges):
        raise NotImplementedError(
            'Use PiecewiseAgeRate or PiecewiseAgeNumber instead of base class PiecewiseAgeParameter')

    def full(self, fullAges, weights = None):
        'Extrapolate from piecewise definition onto an age mesh'

        T = self.makeTransMat(fullAges)
	

        if weights != None:
            # Normalize sum of each column to 1
            
            T /= numpy.outer(numpy.ones(len(fullAges)),
                             T.sum(0))
	    print ("step 2"), T
            # Rescale each column by weights
            T *= numpy.outer(numpy.ones(len(fullAges)),
                             weights)

            # Re-normalize sum of each row to 1
            T /= numpy.outer(T.sum(1),
                             numpy.ones(len(self.ages)))
	return numpy.dot(T, self.values)

    def __repr__(self, header = 'Piecewise constant:'):
        # Find maximum lengths of strings
        valueLen = max([len('%g' % v) for v in self.values])
        ageLen   = max([len('%g' % a) for a in self.ages])

        # Set width of strings to maximum length
        valueFmt = '%%-%dg' % valueLen
        ageFmt   = '%%%dg' % ageLen
        
        strArray = [header]
        for i in xrange(len(self.ages)):
            strVal = '\t= ' \
                     + (valueFmt % self.values[i]) \
                     + '\tfor ' \
                     + (ageFmt % self.ages[i]) \
                     + ' <= a'

            if i < len(self.ages) - 1:
                strVal += ' < ' \
                          + (ageFmt % self.ages[i + 1])

            strArray.append(strVal)

        return '\n'.join(strArray)


class PiecewiseAgeNumber(PiecewiseAgeParameter):
    """
    Number with values[0] between ages[0] and ages[1],
                values[1] between ages[1] and ages[2],
                ...
    """
    def __init__(self, values, ages):
        PiecewiseAgeParameter.__init__(self, values, ages)
    
    def makeTransMat(self, fullAges):
        T = numpy.zeros((len(fullAges), len(self.ages)), dtype = float)


        for (i, a) in enumerate(fullAges):
            for (j, b) in enumerate(self.ages):
                if (j + 1 == len(self.ages) or a < self.ages[j + 1]) and \
                        (i + 1 == len(fullAges) or b < fullAges[i + 1]):
                    if i + 1 == len(fullAges) or j + 1 == len(self.ages):
                        T[i, j] = 1.0
                    else:
                        lenSelf = self.ages[j + 1] - b
                        lenFullInSelf = min(fullAges[i + 1],
                                            self.ages[j + 1]) \
                                            - max(a, b)
                        T[i, j] = float(lenFullInSelf) / float(lenSelf)

        assert numpy.all(abs(T.sum(0) - 1.) < 1e-12), \
               (self.ages, fullAges, T)

        return T

    def __repr__(self):
	
        return PiecewiseAgeParameter.__repr__(
            self,
            header = 'Piecewise constant number:')


class PiecewiseAgeRate(PiecewiseAgeParameter):
    """
    Rate with values[0] between ages[0] and ages[1],
              values[1] between ages[1] and ages[2],
              ...
    """
    def __init__(self, values, ages):
        PiecewiseAgeParameter.__init__(self, values, ages)

    def makeTransMat(self, fullAges):
        T = numpy.zeros((len(fullAges), len(self.ages)), dtype = float)
	#go through full ages groups
        for (i, a) in enumerate(fullAges):
	    # go through reduce age groups
            for (j, b) in enumerate(self.ages):
		#(if pos+1 of reduced age group is the entire length 
		# OR the val of fullage < the val of reduced age at pos+1)
		## AND
		#(if pos+1 of full age is the entire length
		## OR val of reduced age < valur of full age at post +1 )
                if (j + 1 == len(self.ages) or a < self.ages[j + 1]) and \
                        (i + 1 == len(fullAges) or b < fullAges[i + 1]):
		    # if the next positin is the entire length of either fullage or self age
                    if i + 1 == len(fullAges) or j + 1 == len(self.ages):
                        T[i, j] = 1.0
                    else:
                        lenFull = fullAges[i + 1] - a
                        lenSelfInFull = min(fullAges[i + 1],
                                            self.ages[j + 1]) \
                                            - max(a, b)
                        T[i, j] = float(lenSelfInFull) / float(lenFull)

        T /= numpy.outer(T.sum(1), numpy.ones(len(self.ages)))
	
        assert numpy.all(abs(T.sum(1) - 1.) < 1e-12)

        return T

    def __repr__(self):
        return PiecewiseAgeParameter.__repr__(
            self,
            header = 'Piecewise constant rate:')



if __name__ == '__main__':
    import unittest

    class PiecewiseAgeTests(unittest.TestCase):
        'Extend with matrix equality test.'
        
        def assertMatrixEqual(self, A, B):
            'Assert that matrices A and B are equal.'
            self.assertTrue(numpy.all(A == B),
                            '%s != %s' % (A, B))

        def assertMatrixAlmostEqual(self, A, B, places = 7):
            'Assert that matrices A and B are equal to some number of decimal places.'
            self.assertMatrixEqual(numpy.round(A, places),
                                   numpy.round(B, places))


    class PiecewiseAgeNumberTestUnweighted(PiecewiseAgeTests):
        nCoarse = PiecewiseAgeNumber([100., 100.], [0, 5])
        nFine = PiecewiseAgeNumber([20., 80., 100.], [0, 1, 5])
        
        def testRefine(self):
            A = self.nCoarse.full(self.nFine.ages)
            B = self.nFine.values
            self.assertMatrixEqual(A, B)

        def testCoarsen(self):
            A = self.nFine.full(self.nCoarse.ages)
            B = self.nCoarse.values 
            self.assertMatrixEqual(A, B)


    class PiecewiseAgeRateTestUnweighted(PiecewiseAgeTests):
        rCoarse = PiecewiseAgeRate([1., 2.], [0, 5])
        rFine1 = PiecewiseAgeRate([1., 1., 2.], [0, 1, 5])
        rFine2 = PiecewiseAgeRate([2, 0.75, 2.], [0, 1, 5])
        rFine3 = PiecewiseAgeRate([0.6, 1.1, 2.], [0, 1, 5])
        
        def testRefine(self):
            A = self.rCoarse.full(self.rFine1.ages)
            B = self.rFine1.values 
            self.assertMatrixEqual(A, B)

        def testCoarsen1(self):
            A = self.rFine1.full(self.rCoarse.ages)
            B = self.rCoarse.values
            self.assertMatrixEqual(A, B)

        def testCoarsen2(self):
            A = self.rFine2.full(self.rCoarse.ages)
            B = self.rCoarse.values
            self.assertMatrixEqual(A, B)

        def testCoarsen3(self):
            A = self.rFine3.full(self.rCoarse.ages)
            B = self.rCoarse.values
            self.assertMatrixEqual(A, B)


    class PiecewiseAgeRateTestWeighted(PiecewiseAgeTests):
        rCoarse = PiecewiseAgeRate([1., 2.], [0, 5])
        rFine = PiecewiseAgeRate([0.75, 2., 2.], [0, 1, 5])
        wFine = [80, 20, 100]
        
        def testCoarsen(self):
            A = self.rFine.full(self.rCoarse.ages,
                                weights = self.wFine)
            B = self.rCoarse.values
            self.assertMatrixEqual(A, B)


    class PiecewiseAgeRateMortalityTest(PiecewiseAgeTests):
        rFine = PiecewiseAgeRate([0.01265, 0.0004 , 0.0006 , 0.00064,
                                  0.00155, 0.00191, 0.00268, 0.00352,
                                  0.00419, 0.0051 , 0.00644, 0.00893,
                                  0.01255, 0.01756, 0.02664, 0.04331,
                                  0.06946, 0.1328 , 0.22903, 0.35669,
                                  0.5019 , 0.63662],
                                 [ 0,  1,  5, 10,
                                  15, 20, 25, 30,
                                  35, 40, 45, 50,
                                  55, 60, 65, 70,
                                  75, 80, 85, 90,
                                  95, 100])
        rCoarseUnweighted = PiecewiseAgeRate([0.01265, 0.0004 , 0.0006 ,
                                              0.00064, 0.00155, 0.00191,
                                              0.00268, 0.00352, 0.00419,
                                              0.0051 , 0.00644, 0.00893,
                                              0.01255, 0.01756, 0.02664,
                                              0.04331, 0.32108333333333333],
                                             [ 0,  1,  5,
                                              10, 15, 20,
                                              25, 30, 35,
                                              40, 45, 50,
                                              55, 60, 65,
                                              70, 75])
        weights = [4328249 * 0.2, 4328249 * 0.8,
                   4348563, 4808718, 5217959, 4864445,
                   5393328, 5750255, 5548908, 5347483,
                   4929694, 4327654, 3555903, 2700802,
                   2079574, 1690083, 1170277,  659657,
                    268061,   83034,   15361,    1492]
        rCoarseWeighted = PiecewiseAgeRate([0.01265, 0.0004 , 0.0006 ,
                                            0.00064, 0.00155, 0.00191,
                                            0.00268, 0.00352, 0.00419,
                                            0.0051 , 0.00644, 0.00893,
                                            0.01255, 0.01756, 0.02664,
                                            0.04331, 0.12219073692309233],
                                           [ 0,  1,  5,
                                            10, 15, 20,
                                            25, 30, 35,
                                            40, 45, 50,
                                            55, 60, 65,
                                            70, 75])
        
        def testCoarsenUnweighted(self):
            A = self.rFine.full(self.rCoarseUnweighted.ages)
            B = self.rCoarseUnweighted.values
            self.assertMatrixAlmostEqual(A, B)

        def testCoarsenWeighted(self):
            A = self.rFine.full(self.rCoarseWeighted.ages,
                                weights = self.weights)
            B = self.rCoarseWeighted.values
            self.assertMatrixAlmostEqual(A, B)

    unittest.main()
