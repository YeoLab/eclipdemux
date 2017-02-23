#!/bin/env python

"""python distribution setting"""

# from distutils.core import setup
from setuptools import Command
from setuptools import find_packages
from setuptools import setup

from os.path import abspath, dirname, join
from subprocess import call


# automatically set the distribution version
# to be the same as the main package version (eclipdemux_package)
# as defined in eclipdemux_package.__ibit.py__
from eclipdemux_package import __version__


# automatically use README.rst file for the distribution documentation
this_dir = abspath(dirname(__file__))
with open(join(this_dir, 'README.md')
          #, encoding='utf-8'
          ) as file:
    long_description = file.read()


# testing support via py.test library
# with coverage reporting
# run the $ python setup.py test
class RunTests(Command):
    """Run all tests."""
    description = 'run tests'
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        """Run all tests!"""
        errno = call(['py.test', '--cov=eclipdemux_package', '--cov-report=term-missing'])
        raise SystemExit(errno)


setup(
    name='eclipdemux_distribution',
    version='0.0.1',


    # packages=['eclipdemux_pkg'],
    packages = find_packages(exclude=['docs', 'extras', 'tests*', '*/test*']),


    # pip install -e .[test]
    extras_require = {
        'test': ['pylint', 'pytest', 'pytest-cov', 'coverage']
    },

    # run from command line, by typing: $ demux
    entry_points = {
        'console_scripts': [
            'demux=eclipdemux_package.demux:main',
        ],
    },

    # TUTORIAL PURPOSES
    # py_modules=['other_module_1', 'eclipdemux_package.other_module_2']
    )
