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



class SpeciesRegistry (object):  # for recording the genotype of different species
    def __init__(self, initial_number_of_species):
        self.species_list = []  # a species list within which there are different species represented by one marker gene and one genotype list
        self.initial_number_of_species = initial_number_of_species
        for species_marker in range(initial_number_of_species):
            self.species_list.append(species_marker)

    def get_species_community(self, microbiome):  # represent the composition of one microbiome through species_marker regardless of the genotype
        return microbiome

