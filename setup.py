#!/usr/bin/env python

from distutils.core import setup

# Get version number from package
exec(open('dnabclib/version.py').read())

setup(
    name='decontam',
    version=__version__,
    description='Decontaminate paired-end DNA sequencing data',
    author='Kyle Bittinger',
    author_email='kylebittinger@gmail.com',
    url='https://github.com/PennChopMicrobiomeProgram',
    packages=['decontamlib'],
    scripts=['scripts/anti_human.py', 'scripts/make_index.py'],
    )
