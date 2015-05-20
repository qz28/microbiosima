#!/usr/bin/env python
import random

import numpy


class SpeciesRegistry (object):  # for recording the genotype of different species
    def __init__(self, initial_number_of_species):
        self.species_list = []  # a species list within which there are different species represented by one marker gene and one genotype list
        self.initial_number_of_species = initial_number_of_species
        for species_marker in range(initial_number_of_species):
            self.species_list.append(species_marker)

    def get_species_community(self, microbiome):  # represent the composition of one microbiome through species_marker regardless of the genotype
        return microbiome


class Selective_SpeciesRegistry (SpeciesRegistry):  # for recording the genotype of different species
    def __init__(self, initial_number_of_species, number_of_genes, number_of_total_genes, fitness_to_bacterial):
        self.species_list = []  # a species list within which there are different species represented by one marker gene and one genotype list
        self.initial_number_of_species = initial_number_of_species
        self.number_of_total_genes = number_of_total_genes
        self.fitness_list = []  # a list recording the fitness for each species itself (not for host) and
        self.fitness_to_bacterial = fitness_to_bacterial  # the fitness of each gene contributing to bacterial
        self.gene_binary_index = [2 ** k for k in range(number_of_total_genes)]  # binary number representing each gene
        for species_marker in range(initial_number_of_species):
            species = []
            # NOTE: This can be a class to avoid index error. Or use hash with meanful keys
            species.append(species_marker)
            gene_list = random.sample(range(number_of_total_genes), number_of_genes)
            species_fitness = 0
            species_genotype = 0
            for gene in gene_list:
                species_fitness += fitness_to_bacterial[gene]
                species_genotype += self.gene_binary_index[gene]
            self.fitness_list.append(1 ** species_fitness)
            species.append(species_genotype)
            self.species_list.append(species)

    def get_gene_pool(self, microbiome):  # get a gene pool from a microbiome (untested)
        gene_pool = numpy.zeros((self.number_of_total_genes,), dtype=numpy.int)
#         gene_pool = numpy.array([0 for k in range(self.number_of_total_genes)])

        for i in range(len(microbiome)):
            if microbiome[i] > 0:
                supplement = self.number_of_total_genes + 2 - len(bin(self.species_list[i][1]))
                gene_binary_number = '0' * supplement + bin(self.species_list[i][1])[2:]
                gene_pool += numpy.array([int(k) for k in gene_binary_number]) * microbiome[i]
        return gene_pool[::-1]

    def get_fitness_selection(self, microbiome):  # this function is used when species acquisition is totally determined by bacterial fitness
        fitness_selection = []
#         index = 0
        for index, abundance in enumerate(microbiome):
            fitness_selection.append(abundance * self.fitness_list[index])
#             index += 1
        total_fitness = float(sum(fitness_selection))
        return [fitness / total_fitness for fitness in fitness_selection]


class HGT_SpeciesRegistry (Selective_SpeciesRegistry):  # for recording the genotype of different species
    def __init__(self, initial_number_of_species, number_of_genes, number_of_total_genes, fitness_to_bacterial):
        super(HGT_SpeciesRegistry, self).__init__(initial_number_of_species, number_of_genes, number_of_total_genes, fitness_to_bacterial)

    def find_species(self, species):  # when one species is created by hgt, this function return the species index
        if species in self.species_list:
            return self.species_list.index(species)
        else:  # a new species will be added to species list first and then return the species index
            new_index = len(self.species_list)
            self.species_list.append(species)
#             gene = 0
            species_fitness = 0
            for gene, index in enumerate(self.gene_binary_index):
                if species[1] | index == species[1]:
                    species_fitness += self.fitness_to_bacterial[gene]  # get the new bacterial fitness
#                 gene += 1
            self.fitness_list = self.fitness_list + [1 ** species_fitness]  # add it to fitness list
            return new_index

    def get_species_community(self, microbiome):  # represent the composition of one microbiome through species_marker regardless of the genotype
        species_community = microbiome[:self.initial_number_of_species]
        for i in range(self.initial_number_of_species, len(microbiome)):
            species_community[self.species_list[i][0]] += microbiome[i]
        return species_community


