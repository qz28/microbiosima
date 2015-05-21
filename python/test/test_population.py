#!/usr/bin/env python
import unittest
from src.species_registry import SpeciesRegistry
from src.population import Population

class TestPopulation(unittest.TestCase):
    def setUp(self):
        species_registry = SpeciesRegistry(5)
        species_registry.species_list = [[1, 7], [2, 11], [3, 13], [4, 14], [5, 19]]
        self.population = Population(species_registry, [0.2] * 5, 3, 5, 0.5, 0.5)


if __name__ == '__main__':
    unittest.main()
