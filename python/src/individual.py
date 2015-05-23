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


import numpy


class Individual(object):
    def __init__(self, environment, number_of_individual_species, species_registry):  # the host
        self.number_of_environmental_species = len(environment)
        self.microbiome = numpy.random.multinomial(number_of_individual_species, environment)  # initial composition is totally determined by initial environment
        self.number_of_microbes_in_host = number_of_individual_species  # how many microbes in each host
        self.species_registry = species_registry

    def __str__(self):
        self.microbiome_sequence = ['0'] * self.number_of_environmental_species
        for i in range(len(self.microbiome)):
            if self.microbiome[i] != 0:
                self.microbiome_sequence[self.species_registry.species_list[i]] = '1'
        return ''.join(self.microbiome_sequence)

