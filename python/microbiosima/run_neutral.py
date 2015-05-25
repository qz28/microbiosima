################################################################################
#
# Copyright (C) 2014, 2015 Qinglong Zeng, Jeet Sukumaran, Steven Wu and Allen Rodrigo
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
################################################################################

import argparse
import sys

from argparser_util import CustomSingleMetavarFormatter
from population import Population
from species_registry import SpeciesRegistry


MICROBIOSIMA_VERSION = "0.8"
MICROBIOSIMA_DESCRIPTION = "Simulates the evolutionary and ecological dynamics of microbiomes within a population of hosts."



def main():

    parser = argparse.ArgumentParser(
        prog="microbiosima", description=MICROBIOSIMA_DESCRIPTION,
        formatter_class=CustomSingleMetavarFormatter)
    version_string = "%(prog)s " + MICROBIOSIMA_VERSION
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
    CONFIG_DEFAULT = [500, 1000, 150, 1000]
    parser.add_argument(
        "-c", "--config",
        nargs=4,
#         metavar=["config", "2", "3", "4"],
        metavar=("Pop", "Micro", "Spec", "Gen"),
        type=int,
        default=CONFIG_DEFAULT,
        help="""Four Parameters in the following orders:
        (1) population size, (2) microbe size,
        (3) number of species, (4) number of generation [default: %(default)s] 
        """
        )

    parser.add_argument(
        "-o", "--obs", metavar='OBS', type=int,
        help="Number generation for observation [default: %(default)s]", default=100)

    parser.add_argument(
        "-r", "--rep", metavar='REP', type=int,
        help="Number of replication [default: %(default)s]", default=1)

    args = parser.parse_args()
    print args
    try:
        if args.config == CONFIG_DEFAULT:
            print "WARNING! defalut config parameters are use."
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


    def run(species_registry, env, pct_env, pct_pool, rep):
        population = Population(species_registry, env, population_size, microbe_size, pct_env, pct_pool)
        prefix = str(rep + 1) + "_"
        sufix = "_E" + str(pct_env) + "_P" + str(pct_pool) + ".txt"
        print "Output 5 result files in the format of: %s[****]%s" % (prefix, sufix)
        file1 = open(prefix + "fixation" + sufix, 'w')
        file2 = open(prefix + "gamma_diversity" + sufix, 'w')
        file3 = open(prefix + "beta_diversity" + sufix, 'w')
        file4 = open(prefix + "sum" + sufix, 'w')
        file6 = open(prefix + "alpha_diversity" + sufix, 'w')

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
