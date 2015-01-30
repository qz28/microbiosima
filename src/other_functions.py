#!/usr/bin/env python
import random
import numpy
import bisect


def different_element(str1,str2):  #calculate how many elements are different in two strings which are of the same length(for the diversity difference calculation)
    diff=0
    for i in range(len(str1)):
        if str1[i]!=str2[i]:
            diff+=1
    return diff
    
def weighted_choice_b(totals):   #weighted selection through binary search (untested)
    rnd=random.random()*totals[-1]
    return bisect.bisect_left(totals, rnd)
    
def addition_of_arrays(x,y,a,b): #for addition of two numpy.arrays which may have different lengths
    try:
        c=a*x+b*y
    except ValueError:
        d=len(a)-len(b)
        if d>0:
            b=numpy.array(b.tolist()+[0]*d)
        else:
            a=numpy.array(a.tolist()+[0]*(-d))
        c=a*x+b*y
    return c
