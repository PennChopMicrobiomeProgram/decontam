#!/usr/bin/env python

from distutils.core import setup

# Get version number from package
exec(open('decontamlib/version.py').read())

setup(
    name='decontam',
    version=__version__,
    description='Decontaminate paired-end DNA sequencing data',
    author='Kyle Bittinger',
    author_email='kylebittinger@gmail.com',
    url='https://github.com/PennChopMicrobiomeProgram',
    packages=['decontamlib'],
    scripts=['scripts/decontaminate.py', 'scripts/make_index.py', 'scripts/preprocess_report.py'],
    install_requires=["pysam", "biopython"],
    )
