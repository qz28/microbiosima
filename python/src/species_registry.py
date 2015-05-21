#!/usr/bin/env python


class SpeciesRegistry (object):  # for recording the genotype of different species
    def __init__(self, initial_number_of_species):
        self.species_list = []  # a species list within which there are different species represented by one marker gene and one genotype list
        self.initial_number_of_species = initial_number_of_species
        for species_marker in range(initial_number_of_species):
            self.species_list.append(species_marker)

    def get_species_community(self, microbiome):  # represent the composition of one microbiome through species_marker regardless of the genotype
        return microbiome

