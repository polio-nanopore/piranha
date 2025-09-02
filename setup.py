from setuptools import setup, find_packages

setup(
    name="piranha",
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
      package_data={"piranha":["data/*"]}
)