Microbiosima is a collection of python scripts that simulates the evolutionary and ecological
dynamics of microbiomes within a population of hosts.

*****************************************************************************************
The simulation project can be divided into several parts, and currently we have         
finished the neutral framework of microbiome evolution.                                 
                                                                                        
Under our neutral folder, the whole system are consisted of species_registry.py,        
individual.py, population three classes.py. Each class were simulated based             
on the biological meaning in natural communities, such as the construction of           
microbial genomes, the generation of offspring host population, the parental            
inheritance and environmental acquistion of microbes.                                   
                                                                                        
Species_registry: initialize the information of all the microbial species in our system 
                  provide records of the genotypes of all the microbial species         
Individual:       simulate a individual host and the microbial community associated     
                  with the host                                                         
Population:       simulate a population of hosts and their collection of microbiomes,   
                  and all the process that involves in alternation of host-associated   
                  communities such as, the substitution by new offspring, the parental  
                  inheritance of microbial communities, the environmental acquistion of 
                  microbiome communities and the environmental composition affected by  
                  hosts.                                                                
run:              the main function to simulate all the processes with the class         
                  mentioned above.                                                      
microbiosima:     combine all the py files into one.                                    
*****************************************************************************************
                                                                                        
Under our neutral model, several parameters are adjustable:                             
number_of_individuals   (the population size of hosts)                                    
number_of_generations   (the number of host generations that will be simulated)           
number_of_microbes      (number of microbes associated with one host)              
number_of_total_species (the total number of species in the environment)                
env_facotr              (the proportion of environmental_acqusition,                                 
                         1-env_factor is the proportion of parental inheritance)                     
pooled_or_fixed         (the proportion of pooled environmental component in the environment) 
number_generation_
for_observation         (every this many generations the diversities and 
                         other summary statistics are calculated, when it is equal to 1,
                         the observation are made every generation.)
replication             (the number of simulation with the same parameters you want to repeat)  
                                                                                        
our selection model are still under developing process.

Use command "python run(neutral).py env_factor, pooled_or_fixed,number_of_individuals,number_of_microbes,
number_of_total_species,number_of_generations,number_generation_for_observation,replication" to run the
simulation from terminal with eight arguments taken.
You can only provide two arguments for env_factor and pooled_or_fixed, the default settings for other
parameters are following:

number_of_individuals=500
number_of_microbes=1000
number_of_total_species=150
number_of_generations=10000
number_generation_for_observation=100
replication=1
                               
                                                                                        
*****************************************************************************************

Our python scripts are written under python 2.7