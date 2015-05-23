################################################################################
#
# Copyright (C) 2014, 2015 Qinglong Zeng, Jeet Sukumaran, Steven Wu and Allen Rodrigo
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

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



