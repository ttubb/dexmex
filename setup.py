#!/usr/bin/env python

from setuptools import setup, find_packages


setup(
    name='dexmex',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
    ],
    scripts=[
        'dexmex/main.py',
    ],
    include_package_data=True,
    package_data={
        'synum': ['*.tsv']
    },
    entry_points={
        'console_scripts': [
            'dexmex=dexmex.main:main'
        ],
    },
    # Metadata
    author='ttubb',
    author_email='t.tubbesing@uni-bielefeld.de',
    description='Computing local differential expression in metatranscriptomes',
    keywords='metatranscriptome, metagenome, MAG, transcriptomics, rna, rnaseq, deseq',
    url='https://github.com/ttubb/dexmex',
)