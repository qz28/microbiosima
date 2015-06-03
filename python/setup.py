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

import os
# from setuptools import setup
from distutils.core import setup
from microbiosima import __version__
# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name="microbiosima",
    version=__version__,
    author="Qinglong Zeng, Jeet Sukumaran, Steven Wu and Allen Rodrigo",
    author_email="qinglong.zeng@duke.edu",
    description="Simulates the evolutionary and ecological dynamics of microbiomes within a population of hosts.",
    license="GPLv3+",
    keywords="microbiomes netural model",
    url="https://github.com/qz28/microbiosima",
    packages=['microbiosima' ],
    # package_dir={'microbiosima'},
    test_suite="microbiosima.test",
    long_description=read('../README.md'),
    classifiers=[
	'Intended Audience :: Science/Research'
        'License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)',
	'Programming Language :: Python :: 2.7',
	'Topic :: Scientific/Engineering :: Bio-Informatics'
            ],

)
