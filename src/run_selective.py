#!/usr/bin/env python
from SpeciesRegistry import SpeciesRegistry, Selective_SpeciesRegistry
from Population import Population, Selective_Population
from sys import argv 

try:
    script,x,y,population_size,microbe_size,num_species,number_of_generation,number_generation_for_observation,replication,number_of_genes,number_of_total_genes=argv
    y=float(y)
    x=float(x)
    population_size=int(population_size)
    microbe_size=int(microbe_size)
    num_species=int(num_species)
    number_of_generation=int(number_of_generation)
    number_generation_for_observation=int(number_generation_for_observation)
    replication=int(replication)
    number_of_genes=int(number_of_genes)
    number_of_total_genes=int(number_of_total_genes)
except ValueError:
    script,x,y=argv
    y=float(y)
    x=float(x)
    population_size=500
    microbe_size=1000
    num_species=150
    number_of_generation=10
    number_generation_for_observation=10
    replication=1
    number_of_genes=5
    number_of_total_genes=10
if x<0 or x>1 or y<0 or y>1:
    raise ValueError("x and y must be a number between 0 and 1")
    
    

def run(species_registry,env,env_factor,pooled_or_fixed,gene_fitness,rep):
    population=Selective_Population(species_registry,env,population_size,microbe_size,gene_fitness,env_factor,pooled_or_fixed)
    file1=open(str(rep+1)+"_fixation_"+str(y)+"_"+str(x)+"_"+".txt",'w')
    file2=open(str(rep+1)+"_gamma_diversity_"+str(y)+"_"+str(x)+"_"+".txt",'w')
    file3=open(str(rep+1)+"_beta_diversity_"+str(y)+"_"+str(x)+"_"+".txt",'w')
    file4=open(str(rep+1)+"_sum_"+str(y)+"_"+str(x)+"_"+".txt",'w')
    file5=open(str(rep+1)+"_fitness_"+str(y)+"_"+str(x)+"_"+".txt",'w')
    file6=open(str(rep+1)+"_alpha_diversity_"+str(y)+"_"+str(x)+"_"+".txt",'w')
    while population.number_of_generation<=number_of_generation:
        population.sum_species()
        if population.number_of_generation%number_generation_for_observation==0:
            print >>file1, population.ratio_of_fixation()
            print >>file4, population
            print >>file5, str(population.average_of_fitness())+'\t'+str(population.variance_of_fitness())
            print >>file6, population.alpha_diversity()
            print >>file2, population.measure_biodiversity()
            population.microbiome_sequence_alignment()
            population.segregating_site()
            print >>file3, str(population.neleotide_diversity_pi())+'\t'+str(population.number_of_segregating_site)
        population.get_next_gen()

environment=[1/float(num_species) for i in range(num_species)]
fitness_to_bacterial=[0,0,0,0,0,0,0,0,0,0]
species_registry=Selective_SpeciesRegistry(num_species,number_of_genes,number_of_total_genes,fitness_to_bacterial)
gene_fitness=[0,0,0,0,0,0,0,0,0,0]
pooled_or_fixed=y
env_factor=1-x
for rep in range(replication):
    run(species_registry,environment,env_factor,pooled_or_fixed,gene_fitness,rep)
