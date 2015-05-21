#!/usr/bin/env python

import numpy


class Individual(object):
    def __init__(self, environment, number_of_individual_species, species_registry):  # the host
        self.number_of_environmental_species = len(environment)
        self.microbiome = numpy.random.multinomial(number_of_individual_species, environment)  # initial composition is totally determined by initial environment
        self.number_of_microbes_in_host = number_of_individual_species  # how many microbes in each host
        self.species_registry = species_registry

    def __str__(self):
        self.microbiome_sequence = ['0'] * self.number_of_environmental_species
        for i in range(len(self.microbiome)):
            if self.microbiome[i] != 0:
                self.microbiome_sequence[self.species_registry.species_list[i]] = '1'
        return ''.join(self.microbiome_sequence)

