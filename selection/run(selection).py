#!/usr/bin/env python
import SpeciesRegistry
import Population

script,input_number=argv
y=int(input_number)%11
x=int(input_number)/11

def run(species_registry,env,host_selection,gene_fitness,rep):
    population=Population(species_registry,env,500,1000,gene_fitness,env_par,fix_pool,host_selection)
    file1=open(str(rep)+"_fixation_"+str(y)+"_"+str(x)+"_"+str(hgt_rate)+".txt",'w')
    file2=open(str(rep)+"_biodiversity_"+str(y)+"_"+str(x)+"_"+str(hgt_rate)+".txt",'w')
    file3=open(str(rep)+"_statistics_"+str(y)+"_"+str(x)+"_"+str(hgt_rate)+".txt",'w')
    file4=open(str(rep)+"_sum_"+str(y)+"_"+str(x)+"_"+str(hgt_rate)+".txt",'w')
    file5=open(str(rep)+"_fitness_"+str(y)+"_"+str(x)+"_"+str(hgt_rate)+".txt",'w')
    file6=open(str(rep)+"_alpha_diversity_"+str(y)+"_"+str(x)+"_"+str(hgt_rate)+".txt",'w')
    while population.number_of_generation<=10000:
        population.sum_species()
        print >>file1, population.ratio_of_fixation()
        print >>file2, population.measure_biodiversity()
        print >>file5, str(population.average_of_fitness())+'\t'+str(population.variance_of_fitness())
        if population.number_of_generation%100==0:
            print >>file4, population
            print >>file6, population.alpha_diversity()
            population.microbiome_sequence_alignment()
            population.segregating_site()
            print >>file3, str(population.neleotide_diversity_pi())+'\t'+str(population.number_of_segregating_site)
        population.get_next_gen()

num_species=150
environment=[1/float(num_species) for i in range(num_species)]
fitness_to_bacterial=[-1,-1,-1,-1,-1,1,1,1,1,1]
species_selection=1+0.1*x
random.seed(10)
species_registry=SpeciesRegistry(num_species,5,10,fitness_to_bacterial,species_selection)
random.seed()
gene_fitness=[1,1,1,1,1,-1,-1,-1,-1,-1]
host_selection=1+0.1*y
env_par=0.5
fix_pool=0.5
for rep in range(10):
    run(species_registry,environment,host_selection,hgt_rate,gene_fitness,rep)
