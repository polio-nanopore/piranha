from setuptools import setup, find_packages
from setuptools.command.build_py import build_py
import glob
import os
import pkg_resources

from piranha import __version__, _program


setup(name='piranha',
      version=__version__,
      packages=find_packages(),
      scripts=[
            "piranha/scripts/piranha_preprocessing.smk",
            "piranha/scripts/piranha_main.smk",
            "piranha/scripts/piranha_haplotype.smk",
            "piranha/scripts/piranha_consensus.smk",
            "piranha/scripts/piranha_curate.smk",
            "piranha/scripts/piranha_variation.smk",
            "piranha/scripts/piranha_phylo.smk"
            ],
      package_data={"piranha":["data/*"]},
      install_requires=[
            "mako==1.2",
            "pandas~=1.5",
            "snipit>=1.2",
            "biopython",
            "medaka~=1.12.0",
            "numpy~=1.26.4",
            "scipy~=1.11"
      ],
      description='piranha: Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis',
      url='https://github.com/polio-nanopore/piranha',
      author='Aine OToole',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = piranha.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
