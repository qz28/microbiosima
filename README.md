# MICROBIOSIMA

A collection of python scripts that simulates the evolutionary and ecological dynamics of microbiomes within a population of hosts.


## Requirement
* [**python**](https://www.python.org/) 2.7
  * This program is not tested with python 3.0+
* [**NumPy**](http://www.numpy.org/) 1.8.2+



##Components of python scripts

The simulation project can be divided into several parts, and currently we have
finished the neutral framework of microbiome evolution.

- SpeciesRegistry
  - Initialize the information of all the microbial species in our system.
  - provides records of the genotypes of all the microbial species


- Individual:
  - simulate a individual host
  - keep records of the microbial abundance information within the host


- Population
  - simulate a population of hosts and their collection of microbiomes
  - simulate all the process that involves in alternation of host-associated communities
    1. the substitution by new offspring
    2. the parental inheritance of microbial communities
    3. the environmental acquistion of microbiome communities
    4. the environmental composition affected by hosts


##Parameters

Under our neutral model, several parameters are adjustable:

  1. pct_evn: Percentage of environmental acquisition, 1-pct_evn is the proportion of parental inheritance.
  2. pct_pool: Percentage of pooled component in the environment.
  3. population size: the population size of hosts.
  4. microbe size: number of microbes associated with one host
  5. number of species: the total number of species in the environment
  6. number of generations: the number of host generations that will be simulated
  7. number generation for observation: every this many generations the diversities and other summary statistics are calculated
  8. replication: the number of simulation with the same parameters you want to repeat  

##Usage

* Help Menu
```
python src/run_neutral.py -h
```

* Case 1:
```
cd microbiosima
python src/run_neutral.py 0.2 0.5 50 200 20 50
#python src/run_neutral.py arg0 arg1 arg2 arg3 arg4 arg5
```
To run the simulation from terminal with six arguments taken.
  - arg0: pct_env, percentage of environmental acquisition
  - arg1: pct_pool, percentage of pooled component in the environment
  - arg2: population size
  - arg3: microbe size
  - arg4; number of species
  - arg5: number of generations
All parameters are required




* Case 2:


Only two arguments provided for percentage of environmental acquisition, and percentage of pooled component in the environment.
```
cd microbiosima
python src/run_neutral.py 0.2 0.5
#python src/run_neutral.py pct_env pct_pool
#Effectively equals
#python src/run_neutral.py 0.2 0.5 500 1000 150 10000
```
the default settings for other parameters are following:
  - population size=500
  - microbe size=1000
  - number of species=150
  - number of generations=10000


* Additional parameters:
  - `--obs` Number generation for observation [default: 100]
  - `--rep`Number of replication [default: 1]
```
cd microbiosima
python src/run_neutral.py 0.2 0.5 50 200 20 50 --obs 10 --rep 3
```

##Output File format

The format of output filename is "a1_text_a2_a3.txt"
- a1: the number of replicated times
- a2: pct_pool, percentage of pooled component in the environment
- a3: pct_env, percentage of environmental acquisition
- text: category of outputted data
  - alpha_diversity: the average diversity within one host
  - gamma_diversity: the overall diversity within the host population
  - beta_diversity: diversity difference among hosts
  - fixation: fixation time (only ofr pure pooled environment and parental inheritance)
  - sum: the relative abundances of microbial species within population

##Develop

Our selection and HGT model are still under developing process.
