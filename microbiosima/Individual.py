#!/usr/bin/env python
import random
import numpy


class Individual:
    def __init__(self,environment,number_of_individual_species,species_registry,gene_fitness): # the host
        self.number_of_environmental_species=len(environment)
        self.microbiome=numpy.random.multinomial(number_of_individual_species,species_registry.get_fitness_selection(environment)) # initial composition is totally determined by initial environment
        self.number_of_individual_species=number_of_individual_species # how many microbes in each host
        self.gene_pool=species_registry.get_gene_pool(self.microbiome)
        self.species_registry=species_registry
        self.fitness=sum(numpy.array(gene_fitness)*self.gene_pool)/float(number_of_individual_species) # fitness contributing to the host
        self.weighted_fitness=2**(self.fitness)
        
    def __str__(self):
        self.microbiome_sequence=['0']*self.number_of_environmental_species
        for i in range(len(self.microbiome)):
            if self.microbiome[i]!=0:
                self.microbiome_sequence[self.species_registry.species_list[i][0]]='1'
        return ''.join(self.microbiome_sequence)

