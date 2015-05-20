#!/usr/bin/env python
import unittest
import numpy
from src.species_registry import SpeciesRegistry

class TestSR(unittest.TestCase):
    def setUp(self):
        self.species_registry=SpeciesRegistry(5,3,5,[1,0,0,0,-1])

    def test_initial(self):
        self.assertEqual(self.species_registry.gene_binary_index,[1,2,4,8,16])
        count_species={}
        count_fitness={}
        for i in range(1000):
            self.species_registry.__initial__()
            for species in self.species_registry.species_list:
                try:
                    count_species[species[1]]+=1
                except KeyError
                    count_species[species[1]]=1
            for fitness in self.species_registry.species_list:
                try:
                    count_fitness[fitness]+=1
                except KeyError
                    count_fitness[fitness]=1                
            for key in count_species:
                assert(count_species[key]>90&&count_species[key]<110)
            assert(count_fitness[2]>280&&count_fitness[2]<320)
            assert(count_fitness[0.5]>280&&count_fitness[0.5]<320)
            assert(count_fitness[1]>370&&count_fitness[1]<430)

    def test_find_species(self):
        self.species_registry.species_list=[[1,7],[2,11],[3,13],[4,14],[5,19]]
        self.assertEqual(self.species_registry.find_species([3,13]),3)
        self.assertEqual(self.species_registry.find_species([4,14]),4)
        self.assertEqual(self.species_registry.find_species([3,14]),6)
        self.assertEqual(self.species_registry.species_list,[[1,7],[2,11],[3,13],[4,14],[5,19],[3,14]])

    def test_get_species_community(self):
        self.species_registry.species_list=[[1,7],[2,11],[3,13],[4,14],[5,19]]
        self.assertEqual(self.species_registry.get_species_community([1,1,1,1,3]),[1,1,1,1,3])
        self.assertEqual(self.species_registry.get_species_community([1,2,1,1,3]),[1,2,1,1,3])
        self.species_registry.find_species([4,13])
        self.species_registry.find_species([4,11])
        self.species_registry.find_species([3,7])
        self.assertEqual(self.species_registry.get_species_community([1,1,1,1,1,1,1,1]),[1,1,2,3,1])
        
    def test_get_gene_pool(self):
        self.species_registry.species_list=[[1,7],[2,11],[3,13],[4,14],[5,19]]
        self.assertEqual(self.species_registry.get_gene_pool([1,1,1,1,1]),[4,4,3,3,1])
        self.assertEqual(self.species_registry.get_gene_pool([1,2,1,1,1]),[5,5,3,4,1])

    def test_get_fitness_selection(self):
        self.species_registry.fitness_list=[2,1,1,0.5,0.5]
        self.assertEqual(self.species_registry.get_fitness_selection([1,0,1,2,1]),[0.5,0,0.25,0.125,0.125])
        self.assertEqual(self.species_registry.get_fitness_selection([0,10,0,1,10]),[0,0.5,0,0.25,0.25])

if __name__=='__main__':
    unittest.main()
