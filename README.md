# MICROBIOSIMA
Noted on 2024/03/21: We recently fixed a small bug affecting our simulation models of microbiome evolution incorporating host and microbial selection. But all the patterns we have published still holds after the bug fix.

Zeng Q, Sukumaran J, Wu S, Rodrigo A (2015) Neutral Models of Microbiome Evolution. PLoS Comput Biol 11(7): e1004365. doi:10.1371/journal.pcbi.1004365

http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004365

Zeng, Q., Wu, S., Sukumaran, J., & Rodrigo, A. (2017). Models of microbiome evolution incorporating host and microbial selection. Microbiome, 5, 1-16.

https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0343-x

Zeng, Q., & Rodrigo, A. (2018). Neutral models of short-term microbiome dynamics with host subpopulation structure and migration limitation. Microbiome, 6, 1-13.

https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-018-0464-x

Simulates the evolutionary and ecological dynamics of microbiomes within a population of hosts.

There are two versions of this program
* [Python version](#python-version)
* [Java version](#java-version)

###Parameters

Under our neutral model, several parameters are adjustable:
  1. pct_evn: Percentage of environmental acquisition, 1-pct_evn is the proportion of parental inheritance.
  2. pct_pool: Percentage of pooled component in the environment.
  3. population size: the population size of hosts.
  4. microbe size: number of microbes associated with one host
  5. number of species: the total number of species in the environment
  6. number of generations: the number of host generations that will be simulated
  7. number generation for observation: every this many generations the diversities and other summary statistics are calculated
  8. replication: the number of simulation with the same parameters you want to repeat  

##Python Version
### Requirement
* [**python**](https://www.python.org/) 2.7
  * This program is not tested with python 3.0+
* [**NumPy**](http://www.numpy.org/) 1.8.2+



###Usage

##### Help Menu
```bash
python src/run_neutral.py -h
```

##### Case 1

Two arguments provided for percentage of environmental acquisition, and percentage of pooled component in the environment.
```bash
cd microbiosima/python
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

##### Case 2

```bash
cd microbiosima/python
python src/run_neutral.py 0.2 0.5 -c 50 200 20 50
#python src/run_neutral.py pct_env pct_pool -c Pop Micro Spec Gen
```
To run the simulation from terminal with six arguments taken.
- pct_env: percentage of environmental acquisition
- pct_pool: percentage of pooled component in the environment
- Pop: population size
- Micro: microbe size
- Spce: number of species
- Gen: number of generations



##### Additional parameters
  - `--obs` Number generation for observation [default: 100]
  - `--rep`Number of replication [default: 1]

```bash
cd microbiosima
python src/run_neutral.py 0.2 0.5 50 200 20 50 --obs 10 --rep 3
```

###Output File format

The format of output filename is "{a1}_text_E{a2}_P{a3}.txt"
- a1: the number of replicated times
- a2: pct_env, percentage of environmental acquisition
- a3: pct_pool, percentage of pooled component in the environment
- text: category of outputted data
  - alpha_diversity: the average diversity within one host
  - gamma_diversity: the overall diversity within the host population
  - beta_diversity: diversity difference among hosts
  - fixation: fixation time (only for pure pooled environment and parental inheritance)
  - sum: the relative abundances of microbial species within population


### Install
This is not necessary
```bash
# If you have admin right
python setup.py install
# OR at users directory
python setup.py install --user
```


###Components of python scripts

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



##Java Version
### Requirement:
   * [**Java 8**](https://www.java.com/)

### Install
#### Option 1:
Download `microbiosima_*.tar.gz` and uncompress it. Run it from terminal with
```bash
./bin/microbiosima
```

#### Option 2:
Clone git repository and the jar file can be found under the `release` folder.
Alternative, you can recompile it with [ant](http://ant.apache.org/).
```bash
cd microbiosima/java
ant release
cd release/microbiosima*
./bin/microbiosima
```


### Usage
##### Help Menu
```bash
./bin/microbiosima -h
```

##### Case 1

Two arguments provided for percentage of environmental acquisition, and percentage of pooled component in the environment.
```bash
./bin/microbiosima 0.2 0.5
#./bin/microbiosima  pct_env pct_pool
#Effectively equals
#./bin/microbiosima  0.2 0.5 500 1000 150 10000
```
the default settings for other parameters are following:
  - population size=500
  - microbe size=1000
  - number of species=150
  - number of generations=10000

##### Case 2

```bash
./bin/microbiosima 0.2 0.5 -c 50 200 20 50
#./bin/microbiosima pctEnv pctPool -c Pop Micro Spec Gen
```
To run the simulation from terminal with six arguments taken.
- pctEnv: percentage of environmental acquisition
- pctPool: percentage of pooled component in the environment
- Pop: population size
- Micro: microbe size
- Spce: number of species
- Gen: number of generations



##### Additional parameters
  - `--obs` Number generation for observation [default: 100]
  - `--rep`Number of replication [default: 1]




##Development

Our selection and horizontal gene transfer (HGT) model are still under developing process.
