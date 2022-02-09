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
            "piranha/scripts/piranha.smk",
            "piranha/scripts/haplotype_consensus.smk",
            "piranha/scripts/snipit.smk",
            "piranha/scripts/haplotype_reader.smk"
            ],
      package_data={"piranha":["data/*"]},
      install_requires=[
            "mako==1.1.6",
            "snipit==1.0.3"
        ],
      description='piranha: Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis',
      url='https://github.com/aineniamh/piranha',
      author='Aine OToole, Rachel Colquhoun & Corey Ansley',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = piranha.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
