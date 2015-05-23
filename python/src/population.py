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

import math
import random

import numpy

from individual import Individual
from other_functions import addition_of_arrays, different_element


class Population(object):
    def __init__(self, species_registry, environment, number_of_individual,
                 number_of_individual_species, environmental_factor, pooled_or_fixed):
        self.number_of_microbes_in_host = number_of_individual_species  # number of microbes in each host
        self.number_of_environmental_species = len(environment)
        self.number_of_host_in_population = number_of_individual  # number of host in the population
        self.number_of_generation = 0
        self.percentage_of_environmental_acquisition = environmental_factor  # percentage of environmental acquisition
        self.percentage_of_pooled_environmental_component = pooled_or_fixed  # percentage of pooled environmental component
        self.environment = numpy.array(environment)
        self.species_registry = species_registry
        self.composition_of_individual = [Individual(
            environment, number_of_individual_species, species_registry) for _ in range(number_of_individual)]
        self.alpha_diversity = 0
        self.beta_diversity = 0
        self.gamma_diversity = 0
        self.number_of_segregating_site = 0
        self.distance = None
        self.beta_diversity_coef = (self.number_of_host_in_population /
                                    (self.number_of_host_in_population - 1) * 2)


    def sum_species(self):
        self.microbiome_sum = numpy.array([])
        for individual in self.composition_of_individual:
            self.microbiome_sum = addition_of_arrays(1, 1, self.microbiome_sum, individual.microbiome)
        self.microbiome_sum = self.microbiome_sum / float(sum(self.microbiome_sum))  # composition of the microbiomes within the population in terms of genome
        self.species_community = self.microbiome_sum  # composition of the microbiomes within the population in terms of species marker


    def get_from_parent_and_environment(self):  # parental inheritance and environmental acquisition
        parental_contribution = []
        for i in range(self.number_of_host_in_population):  # select parent for next generation
            rnd = random.randint(0, self.number_of_environmental_species - 1)
            parental_contribution.append(
                self.composition_of_individual[rnd].microbiome / float(self.number_of_microbes_in_host))
            self.composition_of_individual[rnd].microbiome *= 2
        environmental_contribution = addition_of_arrays(
            self.percentage_of_pooled_environmental_component,
            1 - self.percentage_of_pooled_environmental_component,
            self.microbiome_sum, self.environment)  # mix pooled and fixed env

        for i in range(self.number_of_host_in_population):  # mix environmental and parental contribution
            mixed_contribution = addition_of_arrays(
                1 - self.percentage_of_environmental_acquisition,
                self.percentage_of_environmental_acquisition,
                parental_contribution[i], environmental_contribution)

            self.composition_of_individual[i].microbiome = numpy.random.multinomial(
                self.number_of_microbes_in_host, mixed_contribution / sum(mixed_contribution))

    def get_next_gen(self):
        self.get_from_parent_and_environment()
        self.pre_species_community = self.species_community
        self.number_of_generation += 1

    def get_distance(self):
        if self.number_of_generation != 0:
            distance_square = 0
            for i in range(self.number_of_environmental_species):
                distance_square += (self.species_community[i] - self.pre_species_community[i]) ** 2
            self.distance = distance_square ** 0.5
        else:
            self.distance = None
        return self.distance

    def ratio_of_fixation(self):
        if (self.percentage_of_pooled_environmental_component == 1 or
                self.percentage_of_environmental_acquisition == 0):
            return (self.species_community.tolist().count(0)) / float(self.number_of_environmental_species)

    def zero_species(self):
        return self.species_community.values().count(0)

    def get_alpha_diversity(self):
        alpha_diversity_list = []
        for i in range(self.number_of_host_in_population):
            individual_species_community = self.species_registry.get_species_community(
                self.composition_of_individual[i].microbiome)
            individual_species_community = individual_species_community / float(sum(individual_species_community))
            biodiversity = 0
            for k in individual_species_community:
                if k > 0 and k < 1:
                    biodiversity = biodiversity - k * math.log(k)
            alpha_diversity_list.append(biodiversity)
        self.alpha_diversity = numpy.average(alpha_diversity_list)
        return self.alpha_diversity

    def get_gamma_diversity(self):  # overall diversity gamma-diversity
        self.gamma_diversity = 0
        for k in self.species_community:
            if k > 0 and k < 1:
                self.gamma_diversity = self.gamma_diversity - k * math.log(k)
        return self.gamma_diversity


    def get_beta_diversity(self):
        # FIXME: you can't call this unless you call microbiome_sequence_alignment() first
        sequence_set = list(set(self.alignment))
        self.beta_diversity = 0
        number_sequence = len(sequence_set)
        if number_sequence > 1:
            for i in range(1, number_sequence):
                fre_i = self.alignment.count(sequence_set[i]) / float(self.number_of_host_in_population)
                for j in range(i):
                    fre_j = self.alignment.count(sequence_set[j]) / float(self.number_of_host_in_population)
                    self.beta_diversity += fre_i * fre_j * different_element(
                        sequence_set[i], sequence_set[j]) / self.number_of_environmental_species
                    # NOTE: Can factor out fre_i, it's faster but less clear
#                     temp_div += fre_j * different_element(
#                         sequence_set[i], sequence_set[j]) / self.number_of_environmental_species
#                 self.beta_diversity += (fre_i * temp_div)

#             self.beta_diversity = self.beta_diversity * self.number_of_host_in_population / (self.number_of_host_in_population - 1) * 2
            self.beta_diversity *= self.beta_diversity_coef
#             return self.beta_diversity
#         else:
        return self.beta_diversity

    def segregating_site(self):
        self.number_of_segregating_site = 0
        for i in range(self.number_of_environmental_species):
            a = self.alignment[0][i]
            for j in range(1, self.number_of_host_in_population):
                if self.alignment[j][i] != a:
                    self.number_of_segregating_site += 1
                    break


    def microbiome_sequence_alignment(self):
        self.alignment = []
        for individual in self.composition_of_individual:
            self.alignment.append(individual.__str__())

    def watterson_tajima(self):  # http://en.wikipedia.org/wiki/Tajima%27s_D
        a_1 = 0
        a_2 = 0
        for i in range(1, self.number_of_host_in_population):
            a_1 += 1 / float(i)
            a_2 += 1 / float(i ** 2)
        self.theta = self.number_of_segregating_site / a_1
        K = 0
        N = 0
        b_1 = (self.number_of_host_in_population + 1) / float(3 * self.number_of_host_in_population - 3)
        b_2 = float(2) / 9 * (self.number_of_host_in_population ** 2 + self.number_of_host_in_population + 3) / (self.number_of_host_in_population ** 2 - self.number_of_host_in_population)
        c_1 = b_1 - 1 / a_1
        c_2 = b_2 - (self.number_of_host_in_population + 2) / (a_1 * self.number_of_host_in_population) + a_2 / (a_1 ** 2)
        e_1 = c_1 / a_1
        e_2 = c_2 / (a_1 ** 2 + a_2)
        for i in range(1, self.number_of_host_in_population):
            for j in range(i):
                K += different_element(self.alignment[i], self.alignment[j])
                N += 1
        k = K / float(N)
        if self.theta != 0:
            self.D = (k - self.theta) / (e_1 * self.number_of_segregating_site + e_2 * self.number_of_segregating_site * (self.number_of_segregating_site - 1)) ** 0.5
        else:
            self.D = None

    def __str__(self):
        return '\t'.join([str(k) for k in self.species_community])

