#!/usr/bin/env python
import unittest
import numpy
from microbiosima.microbiosima import SpeciesRegistry
from microbiosima.microbiosima import Population

class TestSR(unittest.TestCase):
    def setUp(self):
        species_registry=SpeciesRegistry(5,3,5,[1,0,0,0,-1])
        species_registry.species_list=[[1,7],[2,11],[3,13],[4,14],[5,19]]
        self.population=Population(species_registry,[0.2]*5,3,5,[0,0,0,0,0],0,0.5,0.5)

    def test_initial(self):
        self.assertEqual(self.individual.number_of_environmental_species,5)
        count_species=numpy.array([0,0,0,0,0])
        for i in range(1000):
            self.individual.__initial__()
            count_species+=self.individual.microbiome
        for species in count_species:
            assert(species[key]>950&&species[key]<1150)
        self.assertEqual(self.individual.fitness,0)
        self.assertEqual(self.individual.weighted_fitness,1)

    def test_str(self):
        self.individual.microbiome=[1,1,1,1,1]
        self.assertEqual(self.individual.__str__(),"11111")
        self.individual.microbiome=[1,2,2,0,0]
        self.assertEqual(self.individual.__str__(),"11100")
        

if __name__=='__main__':
    unittest.main()
