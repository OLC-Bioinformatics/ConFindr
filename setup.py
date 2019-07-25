#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="confindr",
    version="0.7.1",
    packages=find_packages(),
    entry_points={
       'console_scripts': [
            'confindr.py = confindr_src.confindr:main',
            'confindr = confindr_src.confindr:main',
            'confindr_database_setup = confindr_src.database_setup:main',
            'confindr_create_db = confindr_src.create_genus_specific_db:main'
       ],
    },
    author="Andrew Low",
    author_email="andrew.low@canada.ca",
    url="https://github.com/OLC-Bioinformatics/ConFindr",
    install_requires=['biopython',
                      'pysam',
                      'pytest',
                      'numpy',
                      'rauth']
)
