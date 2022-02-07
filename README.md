# piranha

Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis

<img src="./docs/piranha.svg" width="400">


## Installation instructions
-  Clone the piranha repository with `git clone https://github.com/aineniamh/piranha.git && cd piranha`
-  `conda env create -f environment.yml`
-  `conda activate piranha`
-  `pip install . `

<br>
<h2>Dependencies If <strong>Not</strong> Installing Using the Conda Environment</h2>

<p>
<strong>Note</strong>: Note: we recommend using piranha in the conda environment specified in the environment.yml file as per the instructions above. If you can't use conda for some reason, dependency details can be found in the environment.yml file.
</p>

## Check the install worked
Type (in the <strong>piranha</strong> environment):
	`piranha -v`

## Quick usage

`piranha -i <demultiplexed read directory> -b <path/to/barcodes.csv>`

<br>
<br>
