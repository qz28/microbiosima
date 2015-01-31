# MICOBIOSIMA

A collection of python scripts that simulates the evolutionary and ecological dynamics of microbiomes within a population of hosts.

##components of python scripts

The simulation project can be divided into several parts, and currently we have         
finished the neutral framework of microbiome evolution.                                 
                                                                                                                                                                                                                  
-Species_registry
  *initialize the information of all the microbial species in our system 
  *provide records of the genotypes of all the microbial species
         
-Individual       
  *simulate a individual host
  *keep records of the microbial abundance information within the host
                                                         
-Population
  *simulate a population of hosts and their collection of microbiomes   
  *simulate all the process that involves in alternation of host-associated communities
    1.the substitution by new offspring
    2.the parental inheritance of microbial communities
    3.the environmental acquistion of microbiome communities
    4.the environmental composition affected by hosts
                                                               
-run
  the main function where to start the simulation of all the processes mentioned above.                                                      

##Parameters
                                                                                        
Under our neutral model, several parameters are adjustable:                             
  1.number_of_individuals:the population size of hosts                                    
  2.number_of_generations:the number of host generations that will be simulated           
  3.number_of_microbes:number of microbes associated with one host              
  4.number_of_total_species:the total number of species in the environment                
  5.env_facot:the proportion of environmental_acqusition,1-env_factor is the proportion of parental inheritance                     
  6.pooled_or_fixed:the proportion of pooled environmental component in the environment 
  7.number_generation_for_observation:every this many generations the diversities and other summary statistics are calculated
  8.replication:the number of simulation with the same parameters you want to repeat  

##Usage

Case 1:
Use command "python run(neutral).py arg0,...arg7 " to run the simulation from terminal with eight arguments taken.
  arg0:env_factor
  arg1:pooled_or_fixed
  arg2:number_of_individuals
  arg3:number_of_microbes
  arg4;number_of_total_species
  arg5:number_of_generations
  arg6:number_generation_for_observation
  arg7:replication

Case2:
Only two arguments provided for env_factor and pooled_or_fixed, the default settings for other parameters are following:
  number_of_individuals=500
  number_of_microbes=1000
  number_of_total_species=150
  number_of_generations=10000
  number_generation_for_observation=100
  replication=1

##Output

The format of output filename is "a1_text_a2_a3.txt"
  a1: the number of replicated times
  a2: pooled_or_fixed
  a3: 1-env_factor
  text: category of outputted data
    -alpha_diversity: the average diversity within one host
    -gamma_diversity: the overall diversity within the host population
    -beta_diversity: diversity difference among hosts
    -fixation: fixation time (only ofr pure pooled environment and parental inheritance)
    -sum: the relative abundances of microbial species within population

##Develop
                                                                                       
Our selection and HGT model are still under developing process.
                                                                                                     
Our python scripts are written under python 2.7