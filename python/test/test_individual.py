#!/usr/bin/env python
import unittest
import numpy
from src.species_registry import SpeciesRegistry
from src.individual import Individual

class TestSR(unittest.TestCase):
    def setUp(self):
        self.species_registry = SpeciesRegistry(5)
        self.species_registry.species_list = [[1, 7], [2, 11], [3, 13], [4, 14], [5, 19]]
        self.individual = Individual([0.2] * 5, 5, self.species_registry)

    def test_initial(self):
        self.assertEqual(self.individual.number_of_environmental_species, 5)
        count_species = numpy.array([0, 0, 0, 0, 0])
        for _ in range(1000):
            self.individual = Individual([0.2] * 5, 5, self.species_registry)
            count_species += self.individual.microbiome
        for species in count_species:
            assert(species > 950 and species < 1150)

#     def test_str(self):
#         self.individual.microbiome = [1, 1, 1, 1, 1]
#         self.assertEqual(self.individual.__str__(), "11111")
#         self.individual.microbiome = [1, 2, 2, 0, 0]
#         self.assertEqual(self.individual.__str__(), "11100")


if __name__ == '__main__':
    unittest.main()
