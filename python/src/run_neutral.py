#!/usr/bin/env python
import argparse
import sys

from population import Population
from species_registry import SpeciesRegistry
import numpy

MICROBIOSIMA_VERSION = "0.8"

def main():

    parser = argparse.ArgumentParser(prog="microbiosima")
    version_string = "%(prog)s " + "0.8"
    parser.add_argument(
        "-v", "--version",
        action='version',
        version=version_string,
        )

    parser.add_argument(
        "pct_env", type=float,
        help="Percentage of environmental acquisition")
    parser.add_argument(
        "pct_pool", type=float,
        help="Percentage of pooled environmental component")

    parser.add_argument(
        "config", nargs='*', metavar="config", type=int,
        default=[50, 200, 20, 50],
        help="""Four Parameters in the following orders:
        (1) population size, (2) microbe size,
        (3) number of species, (4) number of generation [default: 50 200 20 50]
        """)

    parser.add_argument(
        "--obs", metavar='', type=int,
        help="Number generation for observation [default: %(default)s]", default=100)

    parser.add_argument(
        "--rep", metavar='', type=int,
        help="Number of replication [default: %(default)s]", default=1)

    args = parser.parse_args()

    print numpy.__version__
    try:
        if args.config == [50, 200, 20, 50]:
            print "WARNING! defalut parameters are use."
        population_size, microbe_size, num_species, number_of_generation = args.config
    except ValueError:
        print "ERROR: config requires exactly 4 arguments! It has %d here %s. EXIT!" % (len(args.config), args.config)
        sys.exit(2)

    number_generation_for_observation = args.obs
    replication = args.rep
    pct_env = args.pct_env
    pct_pool = args.pct_pool

    if pct_env < 0 or pct_env > 1:
        print "ERROR: pct_env (Percentage of environmental acquisition) must be between 0 and 1 (pct_env=%f)! EXIT" % pct_env
        sys.exit(2)
    if pct_pool < 0 or pct_pool > 1:
        print "ERROR: pct_pool (Percentage of pooled environmental component must) must be between 0 and 1 (pct_pool=%f)! EXIT" % pct_pool
        sys.exit(2)

    if replication < 1:
        print "ERROR: Number of replication must be >= 1 (rep=%d)! EXIT" % replication
        sys.exit(2)


    def run(species_registry, env, env_factor, pooled_or_fixed, rep):
        population = Population(species_registry, env, population_size, microbe_size, env_factor, pooled_or_fixed)
        prefix = str(rep + 1)
        sufix = str(pct_pool) + "_" + str(pct_env) + "_" + ".txt"
        print "Output 5 result files in the format of: %s_[****]_%s" % (prefix, sufix)
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




    environment = [1 / float(num_species) for _ in range(num_species)]
    species_registry = SpeciesRegistry(num_species)
    for rep in range(replication):
        run(species_registry, environment, pct_env, pct_pool, rep)
        print "Done replicate %d/%d" % (rep + 1, replication)


if __name__ == "__main__":
    main()
