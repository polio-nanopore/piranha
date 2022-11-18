# piranha

Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis

Piranha is a tool developed to help standardise and streamline sequencing of poliovirus. It's been developed by members of the Rambaut group at the University of Edinburgh as part of the [Poliovirus Sequencing Consortium](https://www.protocols.io/workspaces/poliovirus-sequencing-consortium). Piranha runs an end-to-end read-to-report analysis that produces distributable, interactive reports alongside analysed consensus data. By default, piranha will attempt to generate consensus genomes for populations of the sam It produces an overall report, summarising the entire sequencing run, as well as a sample-specific report. Samples with virus of interest, such as VDPV are highlighted and certain quality-control flags can alert the user if there are issues with the run (such as a failed negative or positive control or identical sequences between samples that may be the result of contamination). 


Any issues or feedback about the analysis or report please flag to this repository.

<img src="./docs/piranha.svg" width="400">

## See example report [here](https://polio-nanopore.github.io/piranha/report.html)

## See example data [here](https://github.com/polio-nanopore/piranha/tree/main/piranha/test/pak_run)


## Installing via ARTIFICE GUI
- Download the release package for your machine from the [ARTIFICE respository](https://github.com/CorwinAnsley/artifice/releases/tag/v1.3.3)

## Installation instructions (quick command line reference)

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

## What to put in a barcodes.csv file?
At a minimum, two columns- one for barcode information and one with the name you'd like your sample to be called. This is also where you can flag which samples are negative or positive controls. Barcode and sample names should be unique and sample names shouldn't contain spaces or special characters. This is because the sample name will get incorporated into the output consensus fasta header and spaces in a fasta header disrupt the sequence ID.

### Minimal example:
```
barcode,sample
barcode01,EDI001
barcode02,EDI002
barcode03,EDI003
barcode04,negative
barcode05,positive
```

You can also include additional information in your barcode.csv. If you include a `date` column or an `EPID` column they'll automatically be included in the fasta header output too. Dates should always be in ISO format (YYYY-MM-DD) for metadata best practice. Piranha has a flag `--all-metadata-to-header` that will take any metadata fields in your barcodes.csv and append them to the final output file (separated by a `|` pipe symbol). Be aware any odd characters or spaces in these fields will also get added to the fasta header and can interfere with downstream phylogenetics you might want to run (e.g. `:`,`;`) can cause issues with some tree-building or reading software). 

### Example with extra information:
```
barcode,sample,EPID,date
barcode01,EDI001,EPI111,2022-10-10
barcode02,EDI002,EPI112,2022-09-20
barcode03,EDI003,EPI113,2022-09,21
barcode04,negative,,
barcode05,positive,,
```

## What should my reads/read directory look like?
Piranha is configured to make analysis as straightforward as possible for users running MinION sequencing. Piranha takes the output of guppy directly and looks for the demultiplexed read directory containing the barcodes you've specified in your barcodes.csv. Guppy outputs directories in the structure:


```
fastq_pass/
├── barcode01/
│    ├── FAT75518_pass_barcode01_6d172881_0.fastq.gz
│    ├── FAT75518_pass_barcode01_6d172881_1.fastq.gz
│    └── FAT75518_pass_barcode01_6d172881_2.fastq.gz
├── barcode02/
│    ├── FAT75518_pass_barcode02_6d172881_0.fastq.gz
│    ├── FAT75518_pass_barcode02_6d172881_1.fastq.gz
│    ├── FAT75518_pass_barcode02_6d172881_2.fastq.gz
│    ├── FAT75518_pass_barcode02_6d172881_3.fastq.gz
│    ├── FAT75518_pass_barcode02_6d172881_4.fastq.gz
│    └── FAT75518_pass_barcode02_6d172881_5.fastq.gz
├── barcode03/
│    ├── FAT75518_pass_barcode03_6d172881_0.fastq.gz
│    ├── FAT75518_pass_barcode03_6d172881_1.fastq.gz
│    └── FAT75518_pass_barcode03_6d172881_2.fastq.gz
└── barcode05/
     ├── FAT75518_pass_barcode05_6d172881_0.fastq.gz
     ├── FAT75518_pass_barcode05_6d172881_1.fastq.gz
     └── FAT75518_pass_barcode05_6d172881_2.fastq.gz
```
This is what piranha will look for. Point the software to the directory containing the different barcodeXX sub-directories and it will iterate within these to find the files. Piranha can accept fastq (fq) or fastq.gz files. It will only attempt to analyse the barcodes present in the input csv file. 


<br>
<br>

## Input configuration

Piranha has been preconfigured with defaults specific to the VP1 protocol developed by the [Polio Sequencing Consortium](https://www.protocols.io/workspaces/poliovirus-sequencing-consortium). All command line arguments (full list below) can be configured either as command line flags when running piranha, or as snake case arguments in a yaml config file (which can then be supplied with the `-c` flag).

For example, you can supply a custom references file (piranha has a default one supplied, which you can access [here](https://github.com/polio-nanopore/piranha/blob/main/piranha/data/references.vp1.fasta)) using the `-r` flag, or by pointing to it within the config file.

<strong>Example Case 1: command line argument</strong>

```(piranha) aine$ piranha -i /localdisk/home/data/minion_run_1/fastq_pass -b /localdisk/home/data/minion_run_1/barcodes.csv -r /localdisk/home/custom_reference_file.fasta ```

This command shows an example piranha run, where everything is in the default settings except the reference file. In this case you point piranha to the read directory, the barcodes.csv file and your custom reference file. 


<strong>Example Case 2: config file</strong>

Example config file (called config.yaml in current working directory). Actual example file can be found [here](https://github.com/polio-nanopore/piranha/blob/main/docs/config.yaml).

```
readdir: /localdisk/home/data/minion_run_1/fastq_pass
barcodes_csv: /localdisk/home/data/minion_run_1/barcodes.csv
reference_sequences: /localdisk/home/custom_reference_file.fasta
```

Then to run piranha you can simply run the command below, and all the information in the file will be included in your run.

```(piranha) aine$ piranha -c config.yaml```

## Flag which samples are controls

Piranha allows you to specify which samples are controls (positive or negative). If the sample name is `negative` or `positive` within the barcode csv file, piranha will automatically detect that these are your controls. See minimal example above for format. 

You can overwrite this if you would rather call your controls something else (like `nc`, `my_fave_control` etc) with the flags `-pc,--positive-control` or `-nc,--negative-control`. 

Example:

```piranha -i path/to/fastq_pass -b barcodes.csv -pc Positive1 -nc "my negative control"```

Alternatively you can supply this in a config file with the fields:
```
positive_control: Positive1
negative_control: Negative1
```

But you need to make sure that the fields match within barcodes.csv. Also note that above because I've put spaces in sample names for my command line example negative control, this command will need quotes around the full name or else the terminal won't interpret it as a single field. ALSO, it's in general better to avoid having spaces in sample names because if you get a consensus sequence out of piranha as a fasta file, record ids are defined as the field up to the first space, so you can lose information in downstream analysis software if you're not careful. Best to just avoid spaces (and also special characters like `:`, `;` and `|`) in general when dealing with this kind of data that might have phylogenetics run on it.

Samples flagged as controls will appear in the report at the end in a separate table as well and will be flagged as either passing (row in table coloured green and a tick appears) or not passing (row in rable coloured red and no tick appears). Piranha's behaviour treats negative controls as passing if there are fewer than the configured minimum number of reads in the sample (Default: 50 reads) and positive controls as passing if there is more than the minimum number of reads in the sample for non-polio enterovirus (Default: 50 reads). 

## Configuring analysis options

Piranha is configured by default to work with the nested amplification VP1 protocol developed by the [Polio Sequencing Consortium](https://www.protocols.io/workspaces/poliovirus-sequencing-consortium). If you're testing another protocol or want to modify the default analysis behaviour you can configure thresholds and models with various command line arguments (the snake case of the long-form arguments can be supplied in an input config file, see example [here](https://github.com/polio-nanopore/piranha/blob/main/docs/config.yaml)). 

### ```--medaka-model```
Piranha uses medaka to generate consensus sequences. Medaka is software developed by Oxford nanopore technologies that explicitly deals with the error profile of nanopore data. Unlike Illumina sequencing, the error profile in sequencing reads generated by nanopore technology is not randomly distributed. Areas of low complexity, repetitive regions and particularly homopolymeric runs are lower in accuracy. Software like medaka (instead of traditional variant calling and consensus generation software) uses machine learning methods to compensate for this non-random error profile and for a long time the leading variant calling software for nanopore sequencing has used machine learning methods (e.g. nanopolish and medaka).

By default the medaka model run is `r941_min_high_g360`. You should ensure to run the appropriate model for your data. The format of medaka model names is of:

```{pore}_{device}_{caller variant}_{caller version}```

so the default model in use assumes you've run 
- a R9.4.1 flow cell 
- on a MinION (or also for GridION can leave it set to `min`- only need to change `min` to `prom` if it was a PromethION run)
- in high accuracy mode (if you've run in fast mode this should be changed)
- with Guppy version 3.6.0

To see the available list of medaka models installed you can use the piranha command 

```(piranha) aine$ piranha --medaka-list-models```

and piranha will check which ones you have installed with your version of medaka, print them to screen and exit. Medaka may lag behind the latest models of guppy available, so use the closest model you can to what you're running but also be aware that the exact model for your version of Guppy may not yet exist if you're on the cutting edge of Guppy versions. Conversely, if you are running a very old version there may also not be an exact medaka model for your data. The developers of medaka suggest to run the correct model for [best results](https://github.com/nanoporetech/medaka#models). They also state: 

>```Where a version of Guppy has been used without an exactly corresponding medaka model, the medaka model with the highest version equal to or less than the guppy version should be selected.```

### Read thresholds

There are a number of default thresholds applied, which can be overwritten depending on your purposes.

```
-n,--min-read-length
-x, --max-read-length
```

The default read length range filters accepts reads between 1000 and 1300 nucleotide bases in length.

```
-d, --min-read-depth
-p, --min-read-pcent
```

These parameters set the minimum number of reads hitting a particular reference in the reference file (and the minimum percentage of reads within the sample) that are necessary to create a binned read group and attempt to make a consensus sequence for that particular sample. By default a minimum of 50 reads are necessary to build a consensus sequence and a minimum of 10% of the sample is required to be represented by that particular reference before it will attempt to create a consensus for this. 

### How many reads should I count as a signal? 

We have set the minimum read depth to be 50 reads in order to attempt to make a consensus. Within piranha, we run minimap2 to map reads against the background reference panel (in a similar manner to [RAMPART](https://github.com/artic-network/rampart)). The top hit within the background reference panel is reported, by default showing the "display_name" field. The categories displayed are:

- Sabin1-Related
- Sabin2-Related
- Sabin3-Related
- WPV1
- WPV2
- WPV3
- NonPolioEV
- Unmapped

When sequencing samples at high depth, using mapping software on the raw nanopore reads (which are error prone) can lead to a certain level of noise. Hits above the minimum read depth threshold (Default >50) are highlighted in red in the final report. If the population of reads mapping to a particular reference successfully makes a consensus sequence at the end of the piranha pipeline, this is an indication of a genuine population of reads rather than noise.

Importantly, the reference that is hit within the background references file does not fully indicate what consensus sequence will be generated from the read population. Further phylogenetic analysis should be performed to confirm the identity of the sequence that you get. 


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
  --all-metadata-to-header
                        Parse all fields from input barcode.csv file and include in the output fasta headers. Be aware spaces in metadata will disrupt the
                        record id, so avoid these.

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
  
