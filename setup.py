#!/usr/bin/env python

from setuptools import setup, find_packages

setup(
    name="confindr",
    version="0.5.0",
    packages=find_packages(),
    entry_points={
       'console_scripts': [
            'confindr.py = confindr_src.confindr:main',
            'confindr_database_setup = confindr_src.database_setup:main'
       ],
    },
    author="Andrew Low",
    author_email="andrew.low@canada.ca",
    url="https://github.com/lowandrew/ConFindr",
    install_requires=['biopython',
                      'pysam',
                      'pytest',
                      'numpy',
                      'rauth']
)
