#!/usr/bin/env python

"""
Setup script for ConFindr.
"""

# Standard imports
from pathlib import Path
import ast
import os

# Third party imports
from setuptools import (
    find_packages,
    setup
)

__author__ = 'adamkoziol'

# Find the version without exec and without distutils
version = {}
version_path = os.path.join(Path(__file__).parent, 'confindr_src', 'version.py')
version_src = Path(version_path).read_text(encoding='utf-8')

# Use AST to safely extract __version__ value
module = ast.parse(version_src)
for node in module.body:
    if isinstance(node, ast.Assign):
        for target in node.targets:
            if isinstance(target, ast.Name) and target.id == '__version__':
                version["__version__"] = ast.literal_eval(node.value)

setup(
    name="confindr",
    version=version['__version__'],
    packages=find_packages(),
    entry_points={
       'console_scripts': [
            'confindr.py = confindr_src.confindr:main',
            'confindr = confindr_src.confindr:main',
            'confindr_database_setup = confindr_src.database_setup:main',
            'confindr_create_db = confindr_src.create_genus_specific_db:main'
       ],
    },
    author="Adam Koziol",
    author_email="adam.koziol@inspection.gc.ca",
    url="https://github.com/OLC-Bioinformatics/ConFindr",
    install_requires=[
       'biopython',
       'pysam',
       'pytest',
       'numpy',
       'rauth',
       'scipy'
      ]
)
