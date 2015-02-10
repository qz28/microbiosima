#!/usr/bin/env python
from sys import argv

from Population import Population
from SpeciesRegistry import SpeciesRegistry


try:
    script, x, y, population_size, microbe_size, num_species, number_of_generation, number_generation_for_observation, replication = argv
    y = float(y)
    x = float(x)
    population_size = int(population_size)
    microbe_size = int(microbe_size)
    num_species = int(num_species)
    number_of_generation = int(number_of_generation)
    number_generation_for_observation = int(number_generation_for_observation)
    replication = int(replication)
except ValueError:
    script, x, y = argv
    y = float(y)
    x = float(x)
    population_size = 500
    microbe_size = 1000
    num_species = 150
    number_of_generation = 10000
    number_generation_for_observation = 10
    replication = 1

if x < 0 or x > 1 or y < 0 or y > 1:
    raise ValueError("x and y must be a number between 0 and 1")



def run(species_registry, env, env_factor, pooled_or_fixed, rep):
    population = Population(species_registry, env, population_size, microbe_size, env_factor, pooled_or_fixed)
    file1 = open(str(rep + 1) + "_fixation_" + str(y) + "_" + str(x) + "_" + ".txt", 'w')
    file2 = open(str(rep + 1) + "_gamma_diversity_" + str(y) + "_" + str(x) + "_" + ".txt", 'w')
    file3 = open(str(rep + 1) + "_beta_diversity_" + str(y) + "_" + str(x) + "_" + ".txt", 'w')
    file4 = open(str(rep + 1) + "_sum_" + str(y) + "_" + str(x) + "_" + ".txt", 'w')
    file6 = open(str(rep + 1) + "_alpha_diversity_" + str(y) + "_" + str(x) + "_" + ".txt", 'w')
    while population.number_of_generation <= number_of_generation:
        population.sum_species()
        if population.number_of_generation % number_generation_for_observation == 0:
            print >> file1, population.ratio_of_fixation()
            print >> file4, population
            print >> file6, population.alpha_diversity()
            print >> file2, population.measure_biodiversity()
            population.microbiome_sequence_alignment()
            population.segregating_site()
            print >> file3, str(population.neleotide_diversity_pi()) + '\t' + str(population.number_of_segregating_site)
        population.get_next_gen()

environment = [1 / float(num_species) for i in range(num_species)]
species_registry = SpeciesRegistry(num_species)
pooled_or_fixed = y
env_factor = x
for rep in range(replication):
    run(species_registry, environment, env_factor, pooled_or_fixed, rep)
