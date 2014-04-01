#!/usr/bin/env python
import random
import numpy
import Individual
import other_functions


class Population:
    def __init__(self,species_registry,environment,number_of_individual,number_of_individual_species,gene_fitness,hgt_rate,environmental_factor,pooled_or_fixed):
        self.number_of_individual_species=number_of_individual_species #number of microbes in each host
        self.number_of_environmental_species=len(environment)
        self.number_of_individual=number_of_individual # number of host in the population
        self.number_of_generation=0
        self.hgt_rate=hgt_rate
        self.environmental_factor=environmental_factor #percentage of environmental acquisition
        self.pooled_or_fixed=pooled_or_fixed #percentage of pooled environmental component
        self.environment=numpy.array(environment) 
        self.gene_fitness=gene_fitness #fitness of each gene contributing to host
        self.species_registry=species_registry
        self.random_binary_index=[2**k for k in range(species_registry.number_of_total_genes)] 
        self.composition_of_individual=[Individual(environment,number_of_individual_species,species_registry,gene_fitness) for k in range(number_of_individual)] 
        self.fitness_collection=[individual.fitness for individual in self.composition_of_individual] #record the host fitness of each generation
        
    def sum_species(self):
        self.microbiome_sum=numpy.array([])
        for individual in self.composition_of_individual:
            self.microbiome_sum=other_functions.addition_of_arrays(1,1,self.microbiome_sum,individual.microbiome)
        self.microbiome_sum=self.microbiome_sum/float(sum(self.microbiome_sum)) #composition of the microbiomes within the population in terms of genome
        self.species_community=self.species_registry.get_species_community(self.microbiome_sum) #composition of the microbiomes within the population in terms of species marker

    def measure_biodiversity(self): #overall diversity gamma-diversity
        self.biodiversity=0 
        for k in self.species_community:
            if k>0 and k<1:
                self.biodiversity=self.biodiversity-k*math.log(k)
        return self.biodiversity
        
    def get_from_parent_and_environment(self): #parental inheritance and environmental acquisition 
        parental_contribution=[]
        self.fitness_collection=[]
        weighted_fitness_totals=[]
        running_totals=0
        for i in range(self.number_of_individual): #preparing for binary search weighted choice
            running_totals+=self.composition_of_individual[i].weighted_fitness
            weighted_fitness_totals.append(running_totals)
        for i in range(self.number_of_individual): #select parent for next generation
            parental_contribution.append(self.composition_of_individual[other_functions.weighted_choice_b(weighted_fitness_totals)].microbiome/float(self.number_of_individual_species))
        environmental_contribution=other_functions.addition_of_arrays(self.pooled_or_fixed,1-self.pooled_or_fixed,self.microbiome_sum,self.environment) #mix pooled and fixed env
        for i in range(self.number_of_individual): #mix environmental and parental contribution
            mixed_contribution=other_functions.addition_of_arrays(1-self.environmental_factor,self.environmental_factor,parental_contribution[i],environmental_contribution)
            self.composition_of_individual[i].microbiome=numpy.random.multinomial(self.number_of_individual_species,self.species_registry.get_fitness_selection(mixed_contribution))
            if self.hgt_rate==0: # when no hgt occurs
                self.composition_of_individual[i].gene_pool=self.species_registry.get_gene_pool(self.composition_of_individual[i].microbiome)
                new_fitness=sum(numpy.array(self.gene_fitness)*self.composition_of_individual[i].gene_pool)/float(self.number_of_individual_species)
                self.fitness_collection.append(new_fitness)
                self.composition_of_individual[i].weighted_fitness=2**new_fitness
                
    def get_from_gene_pool(self):
        for i in range(self.number_of_individual):
            hgt_events=numpy.random.poisson(self.hgt_rate) #determine how many hgt events happen for each host
            if hgt_events>0:
                x=len(self.species_registry.species_list)-len(self.composition_of_individual[i].microbiome) #keep host microbiome list of the same length as species registry
                if x>0:
                    self.composition_of_individual[i].microbiome=numpy.array(self.composition_of_individual[i].microbiome.tolist()+[0]*x)
                gene_pool=self.species_registry.get_gene_pool(self.composition_of_individual[i].microbiome) #gene pool for hgt
                gene_total=0
                cummulative_gene_pool=[]
                for gene in gene_pool:
                    gene_total+=gene
                    cummulative_gene_pool.append(gene_total) # prepare for binary search weighted choice
                species_for_hgt=numpy.random.multinomial(hgt_events,self.composition_of_individual[i].microbiome/float(self.number_of_individual_species))
                species_index=0
                for j in species_for_hgt: # j means how many hgt happens for each species
                    if j>0:
                        old_genotype=self.species_registry.species_list[species_index][1]
                        species_marker=self.species_registry.species_list[species_index][0]
                        while j>0:
                            new_gene=2**(other_functions.weighted_choice_b(cummulative_gene_pool))
                            if old_genotype|new_gene!=old_genotype:
                                random.shuffle(self.random_binary_index)
                                for removed_gene in self.random_binary_index:
                                    if old_genotype|removed_gene==old_genotype: #randomly remove one gene from the genotype
                                        break
                                new_genotype=(old_genotype^removed_gene)|new_gene
                                new_species=[species_marker,new_genotype]
                                new_species_index=self.species_registry.find_species(new_species)
                                self.composition_of_individual[i].microbiome[species_index]=self.composition_of_individual[i].microbiome[species_index]-1 # old species decrease by 1
                                try:
                                    self.composition_of_individual[i].microbiome[new_species_index]=self.composition_of_individual[i].microbiome[new_species_index]+1 # new species increase by 1
                                except IndexError:
                                    self.composition_of_individual[i].microbiome=numpy.array(self.composition_of_individual[i].microbiome.tolist()+[1])
                                if self.composition_of_individual[i].microbiome[species_index]==0: #prevent the number of microbes becoming negative
                                    break
                            j=j-1
                    species_index+=1
            self.composition_of_individual[i].gene_pool=self.species_registry.get_gene_pool(self.composition_of_individual[i].microbiome)
            new_fitness=sum(numpy.array(self.gene_fitness)*self.composition_of_individual[i].gene_pool)/float(self.number_of_individual_species)
            self.fitness_collection.append(new_fitness)
            self.composition_of_individual[i].weighted_fitness=2**new_fitness
            
    def get_next_gen(self):
        self.get_from_parent_and_environment()
        if self.hgt_rate>0:
            self.get_from_gene_pool()
        self.pre_species_community=self.species_community
        self.number_of_generation+=1
        
    def get_distance(self):
        if self.number_of_generation!=0:
            distance_square=0
            for i in range(self.number_of_environmental_species):
                distance_square+=(self.species_community[i]-self.pre_species_community[i])**2
            self.distance=distance_square**0.5
        else:
            self.distance=None
        return self.distance
        
    def ratio_of_fixation(self):
        if self.pooled_or_fixed==1 or self.environmental_factor==0:
            return (self.species_community.tolist().count(0))/float(self.number_of_environmental_species)
            
    def average_of_fitness(self):
        return numpy.average(self.fitness_collection)
        
    def variance_of_fitness(self):
        return numpy.var(self.fitness_collection)    
   
    def zero_species(self):
        return self.species_community.values().count(0)
            
    def alpha_diversity(self):
        alpha_diversity_list=[]
        for i in range(self.number_of_individual):
            individual_species_community=self.species_registry.get_species_community(self.composition_of_individual[i].microbiome)
            individual_species_community=individual_species_community/float(sum(individual_species_community))
            biodiversity=0
            for k in individual_species_community:
                if k>0 and k<1:
                    biodiversity=biodiversity-k*math.log(k)
            alpha_diversity_list.append(biodiversity)
        return numpy.average(alpha_diversity_list)
            
    def microbiome_sequence_alignment(self):
        self.alignment=[]
        for individual in self.composition_of_individual:
            self.alignment.append(individual.__str__())
    
    def neleotide_diversity_pi(self):
        sequence_set=list(set(self.alignment))
        self.pi=0
        number_sequence=len(sequence_set)
        if number_sequence>1:
            for i in range(1,number_sequence):
                for j in range(i):
                    fre_i=self.alignment.count(sequence_set[i])/float(self.number_of_individual)
                    fre_j=self.alignment.count(sequence_set[j])/float(self.number_of_individual)
                    self.pi+=fre_i*fre_j*other_functions.different_element(sequence_set[i],sequence_set[j])/self.number_of_environmental_species
            self.pi=self.pi*self.number_of_individual/(self.number_of_individual-1)
            return self.pi
        else:
            return self.pi
      
    def segregating_site(self):
        self.number_of_segregating_site=0
        for i in range(self.number_of_environmental_species):
            a=self.alignment[0][i]
            for j in range(1,self.number_of_individual):
                if self.alignment[j][i]!=a:
                    self.number_of_segregating_site+=1
                    break
        
    def watterson_tajima(self):
        a_1=0
        a_2=0
        for i in range(1,self.number_of_individual):
            a_1+=1/float(i)
            a_2+=1/float(i**2)
        self.theta=self.number_of_segregating_site/a_1
        K=0
        N=0
        b_1=(self.number_of_individual+1)/float(3*self.number_of_individual-3)
        b_2=float(2)/9*(self.number_of_individual**2+self.number_of_individual+3)/(self.number_of_individual**2-self.number_of_individual)
        c_1=b_1-1/a_1
        c_2=b_2-(self.number_of_individual+2)/(a_1*self.number_of_individual)+a_2/(a_1**2)
        e_1=c_1/a_1
        e_2=c_2/(a_1**2+a_2)
        for i in range(1,self.number_of_individual):
            for j in range(i):
                K+=other_functions.different_element(self.alignment[i],self.alignment[j])
                N+=1
        k=K/float(N)
        if self.theta!=0:
            self.D=(k-self.theta)/(e_1*self.number_of_segregating_site+e_2*self.number_of_segregating_site*(self.number_of_segregating_site-1))**0.5
        else:
            self.D=None

    def __str__(self):
        return '\t'.join([str(k) for k in self.species_community])


