# piranha

Poliovirus Investigation Resource Automating Nanopore Haplotype Analysis


piranha is a tool in development as part of the [Poliovirus Sequencing Consortium](https://polio-nanopore.github.io/). It runs an automated analysis pipeline for sequencing VP1 regions of the poliovirus or non-polio enterovirus genome (whole genome analysis to be implemented soon) and produces an interactive report alongside the consensus data. 

Any issues or feedback about the analysis or report please flag to this repository.

<img src="./docs/piranha.svg" width="400">

## See example report [here](https://polio-nanopore.github.io/piranha/report.html)

## See example data [here](https://github.com/polio-nanopore/piranha/tree/main/piranha/test/pak_run)

## Installation instructions (for in development, will be put on bioconda later)
-  Clone the piranha repository with `git clone https://github.com/polio-nanopore/piranha.git && cd piranha`
-  `conda install -n base -c conda-forge mamba`
-  `mamba env create -f environment.yml`
-  `conda activate piranha`
-  `pip install . `

## Installing via ARTIFICE GUI
- Download the release package for your machine from the [ARTIFICE respository](https://github.com/CorwinAnsley/artifice/releases/tag/v1.3.3)

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

## Check the install worked
Type (in the <strong>piranha</strong> environment):
	`piranha -v`

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

## What should my reads look like?
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
  
