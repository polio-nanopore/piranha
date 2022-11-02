# piranha

Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis

Piranha is a tool developed to help standardise and streamline sequencing of poliovirus. It's been developed by members of the Rambaut group at the University of Edinburgh as part of the [Poliovirus Sequencing Consortium](https://polio-nanopore.github.io/). Piranha runs an end-to-end read-to-report analysis that produces distributable, interactive reports alongside analysed consensus data. It produces an overall report, summarising the entire sequencing run, as well as a sample-specific report. Samples with virus of interest, such as VDPV are highlighted and certain quality-control flags can alert the user if there are issues with the run (such as a failed negative or positive control or identical sequences between samples that may be the result of contamination). 


Any issues or feedback about the analysis or report please flag to this repository.

<img src="./docs/piranha.svg" width="400">

## See example report [here](https://polio-nanopore.github.io/piranha/report.html)



## Installing via ARTIFICE GUI
- Download the release package for your machine from the [ARTIFICE respository](https://github.com/CorwinAnsley/artifice/releases/tag/v1.3.3)

## Installation instructions (quick command reference)

>You need to have Git, a version of conda (link to Miniconda [here](https://docs.conda.io/en/latest/miniconda.html)) and mamba installed to run the following commands. 

Detailed installation instructions given below, but - in brief - to install with mamba run the following in a terminal:

```
git clone https://github.com/polio-nanopore/piranha.git 
cd piranha
mamba env create -f environment.yml
conda activate piranha
pip install . 
```


## Installation instructions (full walkthrough)
If this is your first time running piranha, you’ll need to clone the GitHub repository and install it. If you're running on a Mac machine (OS X) and have never used command line tools such as Git before, you may need to install them. They can be installed easily with the installer. Full instructions can be found [here](https://mac.install.guide/commandlinetools/4.html). Linux systems will have this installed already.

Double check you have git installed by typing the following into a terminal window:

`git --version`

You should see something like the following printed below:

`git version 2.21.1 (Apple Git-122.3)`

If you don't see this, follow the link to the install instructions above to install git. If you're using a Windows machine, it's possible to use Windows Subsystem for Linux, which should have git pre-installed. 

You will also need to have a version of conda installed. We recommend the latest version of Miniconda, which can be accessed for download and install [here](https://docs.conda.io/en/latest/miniconda.html).

I like to keep my repositories in the same place so they’re easy to find, so below I’m making a directory in my home directory (`~`) called repositories and moving into it with `cd`, which is short for "change directory".

```(base) aine$ cd ~```

```(base) aine$ mkdir repositories```

```(base) aine$ cd repositories```

Now we’re ready to clone the piranha GitHub repository. This creates a local copy of the piranha repository on your computer. It also retains a link back to the original repository, so if any updates are made (such as addition of new features or bug fixes), it's very easy to pull those changes down to your local machine and update your version of piranha. More information on git cloning can be found [here](https://git-scm.com/docs/git-clone).


Clone the piranha repository with:

 ```(base) aine$ git clone https://github.com/polio-nanopore/piranha.git```

 ```(base) aine$ cd piranha```

This directory remains linked to the original GitHub copy, so if you need to update it or get any changes you can do that from this location with the following command:

```(base) aine$ git pull```

We now need to create the piranha environment. Hopefully you have mamba installed, check if you do with the following command:

```(base) aine$ mamba –version```

If you see the message `mamba not found` then you should install mamba with:

```(base) aine$ mamba env create -f environment.yml```

Activate the piranha environment:

```(base) aine$ conda activate piranha```

```(piranha) aine$```

The `(piranha)` in the prompt tells you that the piranha environment is activated.

Now you’ll need run install the piranha python package while you’re in the environment:

```(piranha) aine$ pip install .```

The `.` refers to the current working directory (cwd), which should be the piranha repository. To double check you're in the correct directory, you can type `pwd` (print working directory). 

```(piranha) aine$ pwd```

```/localdisk/home/repositories/piranha```

If you see a path printed like the one above, ending with piranha, you know you're in the correct directory. 

Congratulations! You should now have piranha installed. 

## Check the install worked

```(piranha) aine$ piranha --version```

Should return piranha and the version number installed:

```piranha v1.0```

If no errors have come up (such as messages saying "command not found"), you should now be ready to run piranha!

## Installation issue checklist
Sometimes there can be issues unrelated to the commands you've run. 

- Firstly, check you're in the piranha environment. Your prompt should start with `(piranha)`. If not, activate the piranha environment and check your install again.

Example:

```(piranha) aine$ ``` <- correct

```(base) aine$``` <- incorrect

- A common issue can be related to internet connectivity. If a download has failed because of a break in internet, I suggest running through the commands again and just try again. `mamba` and `conda` can cache the files already downloaded, so often you can make progress even if internet connectivity is an issue.
- Similarly, if you have a laptop with not enough storage space, installation and download can fail. To solve this, try clear some space on your machine and try again.
- If you're still having trouble installing via the command line, piranha can be installed using the ARTIFICE GUI (linked above).


## Updating piranha
- Change directory to the piranha repository: `cd piranha`
- Pull the latest changes from GitHub: `git pull`
- Ensure you're in the piranha environment: `conda activate piranha`
- Install the changes: `pip install .`

<br>
<h2>Dependencies If <strong>Not</strong> Installing Using the Conda Environment</h2>

<p>
<strong>Note</strong>: we recommend using piranha in the conda environment specified in the environment.yml file as per the instructions above. If you can't use conda for some reason, dependency details can be found in the environment.yml file.
</p>


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
                        Input config file in yaml format, all command line arguments can be passed via the config file.
  -i READDIR, --readdir READDIR
                        Path to the directory containing fastq read files
  -b BARCODES_CSV, --barcodes-csv BARCODES_CSV
                        CSV file describing which barcodes were used on which sample
  -r REFERENCE_SEQUENCES, --reference-sequences REFERENCE_SEQUENCES
                        Custom reference sequences file.
  -pc POSITIVE_CONTROL, --positive-control POSITIVE_CONTROL
                        Sample name of positive control. Default: `positive`
  -nc NEGATIVE_CONTROL, --negative-control NEGATIVE_CONTROL
                        Sample name of negative control. Default: `negative`

Analysis options:
  -m ANALYSIS_MODE, --analysis-mode ANALYSIS_MODE
                        Specify analysis mode to run. Options: `vp1`. Default: `vp1`
  --medaka-model MEDAKA_MODEL
                        Medaka model to run analysis using. Default: r941_min_high_g360
  --medaka-list-models  List available medaka models and exit.
  -n MIN_READ_LENGTH, --min-read-length MIN_READ_LENGTH
                        Minimum read length. Default: 1000
  -x MAX_READ_LENGTH, --max-read-length MAX_READ_LENGTH
                        Maximum read length. Default: 1300
  -d MIN_READ_DEPTH, --min-read-depth MIN_READ_DEPTH
                        Minimum read depth required for consensus generation. Default: 50
  -p MIN_READ_PCENT, --min-read-pcent MIN_READ_PCENT
                        Minimum percentage of sample required for consensus generation. Default: 10

Output options:
  -o OUTDIR, --outdir OUTDIR
                        Output directory. Default: `analysis-2022-XX-YY`
  -pub PUBLISHDIR, --publishdir PUBLISHDIR
                        Output publish directory. Default: `analysis-2022-XX-YY`
  -pre OUTPUT_PREFIX, --output-prefix OUTPUT_PREFIX
                        Prefix of output directory & report name: Default: `analysis`
  --datestamp DATESTAMP
                        Append datestamp to directory name when using <-o/--outdir>. Default: <-o/--outdir> without a datestamp
  --overwrite           Overwrite output directory. Default: append an incrementing number if <-o/--outdir> already exists
  -temp TEMPDIR, --tempdir TEMPDIR
                        Specify where you want the temp stuff to go. Default: `$TMPDIR`
  --no-temp             Output all intermediate files. For development/ debugging purposes

Misc options:
  --language LANGUAGE   Report language. Options: English, French. Default: English
  --runname RUNNAME     Run name to appear in report. Default: Nanopore sequencing
  --username USERNAME   Username to appear in report. Default: no user name
  --institute INSTITUTE
                        Institute name to appear in report. Default: no institute name
  -t THREADS, --threads THREADS
                        Number of threads. Default: 1
  --verbose             Print lots of stuff to screen
  -v, --version         show program's version number and exit
  -h, --help
  ```
  
  ## pipeline description
  
  - Gathers all read files for a given barcode together
  - Filters these reads by length (Default only reads between 1000 and 1300 bases are included for further analysis)
  - The reads are mapped against a panel of references. By default there is a reference panel included as part of piranha. It includes VP1 sequences from various wild-type polio viruses, reference Sabin-1, Sabin-2 and Sabin-3 sequences for identification of Sabin-like or VDPV polioviruses, and also a number of non-polio enterovirus reference seqeuences. A custom VP1 fasta file can be supplied with the `-r` flag.
  - The read map files are parsed to assign each read to the closest virus sequence in the reference panel. This assigns each read a broad category of either Sabin-1 like, Sabin-2 like, Sabin-3 like, wild-type Poliovirus (1, 2, or 3), non-polio enterovirus, or unmapped.
  - These broad category assignments are used to bin reads for further downstream analysis. Any bin with greater than the minimum read threshold (Default 50 reads, but can be customised) and minimum read percentage (default 10% of sample, but can be customised) is written out in a separate fastq file which will be used to generate the broad-category consensus sequence.
  - For each bin, a consensus sequence is generated using medaka and variation information is calculated for each site in the alignment against the reference. This calculates the consensus variants within each sample. 
  - The variants that are flagged by medaka are assessed for read co-occurance to tease apart variant haplotypes within the sample.
  - For the entire run, and for each individual barcode/ sample, an interactive html report is generated summarising the information. 
  
