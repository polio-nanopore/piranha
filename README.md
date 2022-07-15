# piranha

Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis

<img src="./docs/piranha.svg" width="400">

## See example report [here](https://aineotoole.co.uk/piranha/report.html)

## Installation instructions (for in development, will be put on bioconda later)
-  Clone the piranha repository with `git clone https://github.com/aineniamh/piranha.git && cd piranha`
-  `conda install -n base -c conda-forge mamba`
-  `mamba env create -f environment.yml`
-  `conda activate piranha`
-  `pip install . `

<br>
<h2>Dependencies If <strong>Not</strong> Installing Using the Conda Environment</h2>

<p>
<strong>Note</strong>: we recommend using piranha in the conda environment specified in the environment.yml file as per the instructions above. If you can't use conda for some reason, dependency details can be found in the environment.yml file.
</p>

## Check the install worked
Type (in the <strong>piranha</strong> environment):
	`piranha -v`

## Quick usage

`piranha -i <demultiplexed read directory> -b <path/to/barcodes.csv>`

<br>
<br>

## Full usage
```
usage: 
	piranha -c <config.yaml> [options]
	piranha -i input.csv [options]

Input options:
  -c CONFIG, --config CONFIG
                        Input config file in yaml format, all command line
                        arguments can be passed via the config file.
  -i READDIR, --readdir READDIR
                        Path to the directory containing fastq read files
  -b BARCODES_CSV, --barcodes-csv BARCODES_CSV
                        CSV file describing which barcodes were used on which
                        sample
  -r REFERENCE_SEQUENCES, --reference-sequences REFERENCE_SEQUENCES
                        Custom reference sequences file.

Analysis options:
  -n MIN_READ_LENGTH, --min-read-length MIN_READ_LENGTH
                        Minimum read length.
  -x MAX_READ_LENGTH, --max-read-length MAX_READ_LENGTH
                        Maximum read length.
  -d MIN_READ_DEPTH, --min-read-depth MIN_READ_DEPTH
                        Minimum read depth required for consensus generation.
  -p MIN_READ_PCENT, --min-read-pcent MIN_READ_PCENT
                        Minimum percentage of sample required for consensus
                        generation.

Output options:
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: `analysis-2021-XX-YY`
  -pub PUBLISHDIR, --publishdir PUBLISHDIR
                        Output publish directory. Default: `analysis-2021-XX-
                        YY`
  -pre OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix of output directory & report name: Default:
                        `analysis`
  --datestamp DATESTAMP
                        Append datestamp to directory name when using
                        <-o/--outdir>. Default: <-o/--outdir> without a
                        datestamp
  --overwrite           Overwrite output directory. Default: append an
                        incrementing number if <-o/--outdir> already exists
  -temp TEMPDIR, --tempdir TEMPDIR
                        Specify where you want the temp stuff to go. Default:
                        `$TMPDIR`
  --no-temp             Output all intermediate files. For development/
                        debugging purposes

Misc options:
  -t THREADS, --threads THREADS
                        Number of threads
  --verbose             Print lots of stuff to screen
  -v, --version         show program's version number and exit
  -h, --help
  ```
