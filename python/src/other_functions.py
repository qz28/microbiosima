#!/usr/bin/env python
import bisect
import random

import numpy


def different_element(str1, str2):  # calculate how many elements are different in two strings which are of the same length(for the diversity difference calculation)
    diff = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            diff += 1
    return diff


def weighted_choice_b(totals):  # weighted selection through binary search (untested)
    rnd = random.random() * totals[-1]
    return bisect.bisect_left(totals, rnd)


def addition_of_arrays(weight_a, weight_b, array_a, array_b):  # for addition of two numpy.arrays which may have different lengths
    try:
        array_a *= weight_a
        array_b *= weight_b
        c = array_a * weight_a + array_b * weight_b
    except ValueError:
        d = len(array_a) - len(array_b)
        if d > 0:
            array_b = numpy.array(array_b.tolist() + [0] * d)
        else:
            array_a = numpy.array(array_a.tolist() + [0] * (-d))
        c = array_a * weight_a + array_b * weight_b
    return c



