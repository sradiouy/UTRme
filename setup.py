# -*- coding: utf-8 -*-import re
import os
import sys
import subprocess
try:
    from setuptools import setup, find_packages
except ImportError:
    sys.exit("error: install setuptools")

if sys.version_info < (3, 5):
    raise Exception('utrme requires Python 3.5 or higher.')

def is_tool(name,message):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == os.errno.ENOENT:
            sys.exit(message)
    return True

is_tool("bowtie2","Error: please install bowtie2. See installation instructions.")
is_tool("cutadapt","Error: please install cutadapt. See installation instructions.")
is_tool("bedtools","Error: please install bedtools. See installation instructions.")
is_tool("samtools","Error: please install samtools. See installation instructions.")

setup(name='UTRme',
      version='0.0.3',
      description='UTRme',
      long_description = "UTRme: a scoring based tool to annotate UTR regions in trypanosomatid genomes",
      author='Santiago Radío',
      author_email='sradio91@gmail.com',
      license='MIT',
      packages=['UTRme'],
      entry_points={'console_scripts': ['utrme = UTRme.__main__:main']},
      install_requires=['numpy','xlsxwriter','scipy','python-levenshtein','fuzzywuzzy','regex','pyfaidx','pysam','pandas','BioPython','argparse','twine','matplotlib','seaborn'],
      keywords=['utrme', 'splice-leader', 'polyA','UTR','trans-splicing'],
      classifiers=[
		"Development Status :: 5 - Production/Stable",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		'Operating System :: Unix',
		"Natural Language :: English",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	],
      zip_safe=False
      )
