#!/usr/bin/env python
import unittest
from hgt_abundance 1 import different_element

class TestDE(unittest.TestCase):
    def test(self):
        self.assertEqual(different_element('abs','jsb'),3)
        self.assertEqual(different_element('abs','asb'),2)
        self.assertEqual(different_element('abs','jbs'),1)

if __name__=='__main__':
    unittest.main()