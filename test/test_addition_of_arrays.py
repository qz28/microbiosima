#!/usr/bin/env python
import unittest

import numpy

from src.other_functions import weighted_choice_b


class TestAA(unittest.TestCase):
    def test(self):
        count = [0, 0, 0]
        for i in range(1000):
            count[weighted_choice_b([0.5, 0.7, 1.0])] += 1
        assert count[0] > 450 and count[0] < 550
        assert count[1] > 180 and count[1] < 220
        assert count[2] > 270 and count[2] < 330

if __name__ == '__main__':
    unittest.main()
