#!/usr/bin/env python
import unittest
import numpy
from microbiosima.microbiosima import addition_of_arrays

class TestAA(unittest.TestCase):
    def test(self):
        numpy.array_equal(addition_of_arrays(1,1,numpy.array([1,2,3]),numpy.array([1,2,3,4])),numpy.array([2,4,6,4]))
        numpy.array_equal(addition_of_arrays(1,1,numpy.array([1,1,1]),numpy.array([0,1])),numpy.array([1,2,1]))
        numpy.array_equal(addition_of_arrays(1,1,numpy.array([1,2,3]),numpy.array([1,1,1])),numpy.array([2,3,4]))

if __name__=='__main__':
    unittest.main()
