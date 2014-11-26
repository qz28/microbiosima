Microbiosima is a collection of python scripts that simulates the evolutionary and ecological
dynamics of microbiomes within a population of hosts.

*****************************************************************************************
The simulation project can be divided into several parts, and currently we have         *
finished the neutral framework of microbiome evolution.                                 *
                                                                                        *
Under our neutral model, the whole system are consisted of species_registry,            *
individual_host, population_host three classes. Each class were simulated based         *
on the biological meaning in natural communities, such as the construction of           *
microbial genomes, the generation of offspring host population, the parental            *
inheritance and environmental acquistion of microbes.                                   *
                                                                                        *
Under our neutral model, several parameters are adjustable:                             *
number_of_individuals (the population size of hosts)                                    *
number_of_generations (the number of host generations that will be simulated)           *
number_of_individual_species (number of microbes associated with one host)              *
number_of_total_species (the total number of species in the environment)                *
env_facotr (the proportion of environmental_acqusition,                                 *
            1-env_factor is the proportion of parental inheritance)                     *
pooled_or_fixed (the proportion of pooled environmental component in the environment)   *
                                                                                        *
our selection model are still under developing process                                  *
                                                                                        *
*****************************************************************************************

Our python scripts are written under python 2.7