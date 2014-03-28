#!/usr/bin/env python
import unittest
import numpy
from hgt_abundance import SpeciesRegistry

class TestSR(unittest.TestCase):
    def setUp(self):
        self.species_registry=SpeciesRegistry(5,1,3,[1,0,-1])
        
    def test_initial(self):
        self.assertEqual(self.species_registry.gene_binary_index,[1,2,4])
        for i in range(self.species_registry.initial_number_of_species):
            self.assertEqual(i,self.species_registry.species_list[i][0])
            assert(self.species_registry.species_list[i][1] in self.species_registry.gene_binary_index)
            assert(self.species_registry.fitness_list[i] in [2,1,0.5])
            
    def test_find_species(self):
        self.assertEqual(self.species_registry.find_species(self.species_registry.species_list[3]),3)
        self.assertEqual(self.species_registry.find_species(self.species_registry.species_list[4]),4)
        assert(self.species_registry.find_species([4,1]) in [4,5,6])
        assert(self.species_registry.find_species([4,2]) in [4,5,6])
        assert(self.species_registry.find_species([4,4]) in [4,5,6])
        self.assertEqual(len(self.species_registry.species_list),7)
        
    def test_get_species_community(self):
        self.assertEqual(self.species_registry.get_species_community([1,1,1,1,3]),[1,1,1,1,3])
        self.assertEqual(self.species_registry.get_species_community([1,2,1,1,3]),[1,2,1,1,3])
        self.species_registry.find_species([4,1])
        self.species_registry.find_species([4,2])
        self.species_registry.find_species([4,4])
        self.assertEqual(self.species_registry.get_species_community([1,1,1,1,1,1,1]),[1,1,1,1,3])
        
    def test_get_fitness_selection(self):
        self.assertEqual(self.species_registry.get_fitness_selection([1,0,0,0,0]),[1,0,0,0,0])
        self.species_registry.find_species([4,1])
        self.species_registry.find_species([4,2])
        self.species_registry.find_species([4,4])
        self.assertEqual(self.species_registry.get_fitness_selection([1,0,0,0,0]),[1,0,0,0,0])        

if __name__=='__main__':
    unittest.main()