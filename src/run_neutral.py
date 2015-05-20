#!/usr/bin/env python
from sys import argv

from src import population
from src import species_registry


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
#     print argv
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
    population = population(species_registry, env, population_size, microbe_size, env_factor, pooled_or_fixed)
    prefix = str(rep + 1)
    sufix = str(y) + "_" + str(x) + "_" + ".txt"
    file1 = open(prefix + "_fixation_" + sufix, 'w')
    file2 = open(prefix + "_gamma_diversity_" + sufix, 'w')
    file3 = open(prefix + "_beta_diversity_" + sufix, 'w')
    file4 = open(prefix + "_sum_" + sufix, 'w')
    file6 = open(prefix + "_alpha_diversity_" + sufix, 'w')
    while population.number_of_generation <= number_of_generation:
        population.sum_species()
        if population.number_of_generation % number_generation_for_observation == 0:
            # NOTE: separate summarise and get stats
            # population.calcalute_all_stats()
            # then use accessor to get values.
            print >> file1, population.ratio_of_fixation()
            print >> file4, population
            print >> file6, population.get_alpha_diversity()
            print >> file2, population.get_gamma_diversity()
            population.microbiome_sequence_alignment()
            population.segregating_site()
            print >> file3, str(population.get_beta_diversity()) + '\t' + str(population.number_of_segregating_site)
        population.get_next_gen()


environment = [1 / float(num_species) for i in range(num_species)]
species_registry = species_registry(num_species)
pooled_or_fixed = y
env_factor = x
for rep in range(replication):
    run(species_registry, environment, env_factor, pooled_or_fixed, rep)



"""
parser = argparse.ArgumentParser()
parser.add_argument("-p", help="population_size")
parser.add_argument("-m", help="microbe_size")
parser.add_argument("-g", help="number of generation")
parser.add_argument("-r", help="replication")
parser.add_argument("-d", "--debug", action="count", default=0,
                        help="increase debugging level")
"""
