#!/usr/bin/env python
import unittest
from microbiosima.microbiosima import different_element

class TestDE(unittest.TestCase):
    def test(self):
        self.assertEqual(different_element('abc','abs'),1)
        self.assertEqual(different_element('aaa','bbb'),3)
        self.assertEqual(different_element('xysh','jshs'),4)

if __name__=='__main__':
    unittest.main()
