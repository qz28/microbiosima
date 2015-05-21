#!/usr/bin/env python
import unittest
from src.species_registry import SpeciesRegistry

class TestSR(unittest.TestCase):
    def setUp(self):
        self.species_registry = SpeciesRegistry(5)

    def test_initial(self):

        count_species = {}
        count_fitness = {}
        for i in range(1000):
            self.species_registry = SpeciesRegistry(5)
            for species in self.species_registry.species_list:
                try:
                    count_species[species] += 1
                except KeyError:
                    count_species[species] = 1
            for fitness in self.species_registry.species_list:
                try:
                    count_fitness[fitness] += 1
                except KeyError:
                    count_fitness[fitness] = 1
            for key in count_species:
                expected = i + 1
                self.assertEqual(count_species[key], expected)
                self.assertEqual(count_fitness[key], expected)



    def test_get_species_community(self):
        self.species_registry.species_list = [[1, 7], [2, 11], [3, 13], [4, 14], [5, 19]]
        self.assertEqual(self.species_registry.get_species_community([1, 1, 1, 1, 3]), [1, 1, 1, 1, 3])
        self.assertEqual(self.species_registry.get_species_community([1, 2, 1, 1, 3]), [1, 2, 1, 1, 3])


if __name__ == '__main__':
    unittest.main()
