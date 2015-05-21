#!/usr/bin/env python
import unittest
import numpy
from src.other_functions import weighted_choice_b, different_element, addition_of_arrays


class TestOtherFunctions(unittest.TestCase):
    def test_weighted_choice(self):
        count = [0, 0, 0]
        for _ in range(1000):
            count[weighted_choice_b([0.5, 0.7, 1.0])] += 1
        assert count[0] > 450 and count[0] < 550
        assert count[1] > 180 and count[1] < 220
        assert count[2] > 270 and count[2] < 330

    def test_different_element(self):
        self.assertEqual(different_element('abc', 'abs'), 1)
        self.assertEqual(different_element('aaa', 'bbb'), 3)
        self.assertEqual(different_element('xysh', 'jshs'), 4)

    def test_addation_of_arrays(self):
        numpy.array_equal(addition_of_arrays(1, 1, numpy.array([1, 2, 3]), numpy.array([1, 2, 3, 4])), numpy.array([2, 4, 6, 4]))
        numpy.array_equal(addition_of_arrays(1, 1, numpy.array([1, 1, 1]), numpy.array([0, 1])), numpy.array([1, 2, 1]))
        numpy.array_equal(addition_of_arrays(1, 1, numpy.array([1, 2, 3]), numpy.array([1, 1, 1])), numpy.array([2, 3, 4]))



if __name__ == '__main__':
    unittest.main()
