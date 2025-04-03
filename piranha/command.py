#!/usr/bin/env python3
from piranha import __version__

from piranha.input_parsing import initialising as init

from piranha.utils import misc
from piranha.utils import dependency_checks
from piranha.utils import data_install_checks
from piranha.input_parsing import analysis_arg_parsing
from piranha.input_parsing import directory_setup
from piranha.input_parsing import input_qc
from piranha.analysis import phylo_functions

from piranha.report.make_report import make_output_report
from piranha.utils.log_colours import green,cyan,red
from piranha.utils.config import *

import os
import sys
import shutil
import yaml
import argparse

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):
 
    parser = argparse.ArgumentParser(add_help=False,
    description=misc.preamble(__version__),
    usage='''
\tpiranha -c <config.yaml> [options]
\tpiranha -b <barcodes.csv> -i <demultiplexed fastq_dir> [options]
''')

    i_group = parser.add_argument_group('Input options')
    i_group.add_argument('-c',"--config", action="store",help="Input config file in yaml format, all command line arguments can be passed via the config file.", dest="config")
    i_group.add_argument('-i','--readdir',help="Path to the directory containing fastq read files",dest="readdir")
    i_group.add_argument('-b','--barcodes-csv',help="CSV file describing which barcodes were used on which sample",dest="barcodes_csv")
    i_group.add_argument("-r","--reference-sequences",action="store",dest="reference_sequences",help="Custom reference sequences file.")
    i_group.add_argument("-rg","--reference-group-field",action="store",help=f"Specify reference description field to group references by. Default: `{VALUE_REFERENCE_GROUP_FIELD}`")
    i_group.add_argument("-nc","--negative-control",action="store",help=f"Sample name of negative control. If multiple samples, supply as comma-separated string of sample names. E.g. `sample01,sample02` Default: `{VALUE_NEGATIVE}`")
    i_group.add_argument("-pc","--positive-control",action="store",help=f"Sample name of positive control. If multiple samples, supply as comma-separated string of sample names. E.g. `sample01,sample02`. Default: `{VALUE_POSITIVE}`")
    i_group.add_argument("-pr","--positive-references",action="store",help=f"Comma separated string of sequences in the reference file to class as positive control sequences.")

    analysis_group = parser.add_argument_group('Analysis options')
    analysis_group.add_argument("-s","--sample-type",action="store",help=f"Specify sample type. Options: `stool`, `environmental`, `isolate`. Default: `{VALUE_SAMPLE_TYPE}`")
    analysis_group.add_argument("-m","--analysis-mode",action="store",help=f"Specify analysis mode to run, for preconfigured defaults. Options: `vp1`, `wg`. Default: `{VALUE_ANALYSIS_MODE}`")
    analysis_group.add_argument("--medaka-model",action="store",help=f"Medaka model to run analysis using. Default: {VALUE_DEFAULT_MEDAKA_MODEL}")
    analysis_group.add_argument("--medaka-list-models",action="store_true",help="List available medaka models and exit.")
    analysis_group.add_argument("-q","--min-map-quality",action="store",type=int,help=f"Minimum mapping quality. Range 0 to 60, however 0 can imply a multimapper. Default: {VALUE_MIN_MAP_QUALITY}")
    analysis_group.add_argument("-n","--min-read-length",action="store",type=int,help=f"Minimum read length. Default: {READ_LENGTH_DICT[VALUE_ANALYSIS_MODE][0]}")
    analysis_group.add_argument("-x","--max-read-length",action="store",type=int,help=f"Maximum read length. Default: {READ_LENGTH_DICT[VALUE_ANALYSIS_MODE][1]}")
    analysis_group.add_argument("-d","--min-read-depth",action="store",type=int,help=f"Minimum read depth required for consensus generation. Default: {VALUE_MIN_READS}")
    analysis_group.add_argument("-p","--min-read-pcent",action="store",type=float,help=f"Minimum percentage of sample required for consensus generation. Default: {VALUE_MIN_PCENT}")
    analysis_group.add_argument("-a","--min-aln-block",action="store",type=float,help=f"Minimum alignment block length. Default: 0.6*MIN_READ_LENGTH")
    analysis_group.add_argument("--primer-length",action="store",type=int,help=f"Length of primer sequences to trim off start and end of reads. Default: {VALUE_PRIMER_LENGTH}")
    analysis_group.add_argument("-mo","--minimap2-options",action="store",nargs='*',
                                help=f"Specify any number of minimap2 options to overwrite the default mapping settings. The options take the form `flag=value` and can be any number of space-delimited options. Default: {VALUE_DEFAULT_MINIMAP2}. Example: For short reads of a sample diverged from the reference, we suggest using `-mo k=5 w=4`, which overwrites the minimap2 option `-x asm20`.")

    haplo_group = parser.add_argument_group('Haplotyping options')
    haplo_group.add_argument("-rh","--run-haplotyping",action="store_true",help=f"Trigger the optional haplotyping module. Additional dependencies may need to be installed.")
    haplo_group.add_argument("-hs","--haplotype-sample-size",action="store",type=int,help=f"Number of reads to downsample to for haplotype calling. Default: {VALUE_HAPLOTYPE_SAMPLE_SIZE}")
    haplo_group.add_argument("-hf","--min-allele-frequency",action="store",type=float,help=f"Minimum allele frequency to call. Note: setting this below 0.07 may significantly increase run time. Default: {VALUE_MIN_ALLELE_FREQUENCY}")
    haplo_group.add_argument("-hx","--max-haplotypes",action="store",type=int,help=f"Maximum number of haplotypes callable within reference group. Default: {VALUE_MAX_HAPLOTYPES}")
    haplo_group.add_argument("-hdist","--min-haplotype-distance",action="store",type=int,help=f"Minimum number of SNPs between haplotypes. Default: {VALUE_MIN_HAPLOTYPE_DISTANCE}")
    haplo_group.add_argument("-hd","--min-haplotype-depth",action="store",type=int,help=f"Minimum number of reads in a given haplotype. Default: {VALUE_MIN_HAPLOTYPE_DEPTH}")

    phylo_group = parser.add_argument_group('Phylogenetics options')
    phylo_group.add_argument("-rp","--run-phylo",action="store_true",help=f"Trigger the optional phylogenetics module. Additional dependencies may need to be installed.")
    phylo_group.add_argument("-sd","--supplementary-datadir",action="store",help=f"Path to directory containing supplementary sequence FASTA file and optional metadata to be incorporated into phylogenetic analysis.")
    phylo_group.add_argument("-pcol","--phylo-metadata-columns",action="store",help=f"Columns in the barcodes.csv file to annotate the phylogeny with. Default: {VALUE_PHYLO_METADATA_COLUMNS}")
    phylo_group.add_argument("-smcol","--supplementary-metadata-columns",action="store",help=f"Columns in the supplementary metadata to annotate the phylogeny with. Default: {VALUE_SUPPLEMENTARY_METADATA_COLUMNS}")
    phylo_group.add_argument("-smid","--supplementary-metadata-id-column",action="store",help=f"Column in the supplementary metadata files to match with the supplementary sequences. Default: {VALUE_SUPPLEMENTARY_METADATA_ID_COLUMN}")
    phylo_group.add_argument("-ud","--update-local-database",action="store_true",help=f"Add newly generated consensus sequences (with a distance greater than a threshold (--local-database-threshold) away from Sabin, if Sabin-related) and associated metadata to the supplementary data directory.")
    phylo_group.add_argument("-dt","--local-database-threshold",action="store_true",help=f"The threshold beyond which Sabin-related sequences are added to the supplementary data directory if update local database flag used. Default: {VALUE_LOCAL_DATABASE_THRESHOLD}")

    o_group = parser.add_argument_group('Output options')
    o_group.add_argument('-o','--outdir', action="store",help=f"Output directory. Default: `{VALUE_OUTPUT_PREFIX}-2022-XX-YY`")
    o_group.add_argument('-pub','--publishdir', action="store",help=f"Output publish directory. Default: `{VALUE_OUTPUT_PREFIX}-2022-XX-YY`")
    o_group.add_argument('-pre','--output-prefix',action="store",help=f"Prefix of output directory & report name: Default: `{VALUE_OUTPUT_PREFIX}`",dest="output_prefix")
    o_group.add_argument('--datestamp', action="store",help="Append datestamp to directory name when using <-o/--outdir>. Default: <-o/--outdir> without a datestamp")
    o_group.add_argument('--overwrite', action="store_true",help="Overwrite output directory. Default: append an incrementing number if <-o/--outdir> already exists")
    o_group.add_argument('-temp','--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: `$TMPDIR`")
    o_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files. For development/ debugging purposes",dest="no_temp")
    o_group.add_argument("-ff","--fasta-header-fields",action="store",help=f"Comma-separated string of header fields wanted that can be sourced from input csv file. E.g. `sample,barcode,runname`.")
    o_group.add_argument('--all-metadata-to-header',action="store_true",dest=KEY_ALL_METADATA,help="Parse all fields from input barcode.csv file and include in the output fasta headers. Be aware spaces in metadata will disrupt the record id, so avoid these.")
    o_group.add_argument('--language',action="store",help=f"Output report language. Options: English, French. Default: {VALUE_LANGUAGE}")
    o_group.add_argument('--save-config',action="store_true",help=f"Output the config file with all parameters used")
    o_group.add_argument('--archive-fastq',action="store_true",help=f"Write the supplied fastq_pass directory to the output directory.")
    o_group.add_argument('--archivedir',action="store",help=f"Configure where to put the fastq_pass files, default in the output directory.")

    misc_group = parser.add_argument_group('Misc options')
    misc_group.add_argument('--runname',action="store",help=f"Run name to appear in report. Default: {VALUE_RUNNAME}")
    misc_group.add_argument('--username',action="store",help="Username to appear in report. Default: no user name")
    misc_group.add_argument('--institute',action="store",help="Institute name to appear in report. Default: no institute name")
    misc_group.add_argument('--notes',action="store",help="Miscellaneous notes to appear at top of report. Default: no notes")
    misc_group.add_argument('--orientation',action="store",help="Orientation of barcodes in wells on a 96-well plate. If `well` is supplied as a column in the barcode.csv, this default orientation will be overwritten. Default: `vertical`. Options: `vertical` or `horizontal`.")
    misc_group.add_argument('-t', '--threads', action='store',dest="threads",type=int,help="Number of threads. Default: 1")
    misc_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    misc_group.add_argument("-v","--version", action='version', version=f"piranha {__version__}")
    misc_group.add_argument("-h","--help",action="store_true",dest="help")

    """
    Exit with help menu if no args supplied
    """

    if len(sysargs)<1: 
        parser.print_help()
        sys.exit(0)
    else:
        args = parser.parse_args(sysargs)
        if args.help:
            parser.print_help()
            sys.exit(0)

    dependency_checks.check_dependencies(DEPENDENCY_LIST, MODULE_LIST)

    # Initialise config dict
    config = init.setup_config_dict(cwd,args.config)
    
    # Check if valid sample_type
    analysis_arg_parsing.sample_type(args.sample_type,config)

    # Checks access to package data
    data_install_checks.check_install(args.language,config)

    analysis_arg_parsing.analysis_mode(args.analysis_mode,config)

    # orientation for pcr plate (for report figure)
    misc.add_check_valid_arg(KEY_ORIENTATION,args.orientation,VALID_ORIENTATION,config)

    # grabs the snakefile
    snakefile = data_install_checks.get_snakefile(thisdir,"main")
    
    # Checks medaka options if non default values used.
    analysis_arg_parsing.medaka_options_parsing(args.medaka_model,args.medaka_list_models,config)

    # Configures which analysis options to run
    analysis_arg_parsing.analysis_group_parsing(args.min_read_length,
                                                args.max_read_length,
                                                args.min_read_depth,
                                                args.min_read_pcent,
                                                args.min_aln_block,
                                                args.primer_length,
                                                args.min_map_quality,
                                                config)

    config[KEY_MINIMAP2_OPTIONS] = analysis_arg_parsing.minimap2_options_parsing(args.minimap2_options,config)

    # Configures which haplotype defaults and whether to run haplo calling
    analysis_arg_parsing.haplo_group_parsing(args.run_haplotyping,
                        args.haplotype_sample_size,
                        args.min_allele_frequency,
                        args.max_haplotypes,
                        args.min_haplotype_distance,
                        args.min_haplotype_depth,
                        config)

    misc.add_arg_to_config(KEY_ALL_METADATA,args.all_metadata_to_header,config)

    input_qc.parse_input_group(args.barcodes_csv,
                                args.readdir,
                                args.reference_sequences,
                                args.reference_group_field,
                                config)

    input_qc.control_group_parsing(args.positive_control,
                                    args.negative_control,
                                    args.positive_references,
                                    config)

    # sets up the output dir, temp dir, and data output desination
    directory_setup.output_group_parsing(args.outdir,
                                        args.output_prefix,
                                        args.overwrite,
                                        args.datestamp,
                                        args.tempdir,
                                        args.no_temp,
                                        args.archive_fastq,
                                        args.archivedir,
                                        config)

    init.configure_fasta_header(args.fasta_header_fields,config)

    init.misc_args_to_config(args.verbose,
                                args.threads,
                                args.username,
                                args.institute,
                                args.runname,
                                args.notes,
                                config)
    
    # runs qc checks on the phylo input options and configures the phylo settings
    # now need tempdir for this parsing, so run after directory_setup
    # also needs runname to not add runname.today.fasta to the db
    input_qc.phylo_group_parsing(args.run_phylo, 
                                args.update_local_database,
                                args.supplementary_datadir,
                                args.phylo_metadata_columns,
                                config[KEY_BARCODES_CSV],
                                args.supplementary_metadata_columns,
                                args.supplementary_metadata_id_column,
                                args.local_database_threshold,
                                config)

    if config[KEY_RUN_PHYLO]:
        # checks the phylo-specific dependencies
        dependency_checks.check_dependencies(PHYLO_DEPENDENCY_LIST, PHYLO_MODULE_LIST)

    # ready to run? either verbose snakemake or quiet mode
    init.set_up_verbosity(config)

    preprocessing_snakefile = data_install_checks.get_snakefile(thisdir,"preprocessing")
    phylo_snakefile = data_install_checks.get_snakefile(thisdir,"phylo")
    haplo_snakefile = data_install_checks.get_snakefile(thisdir,"haplotype")

    # output an optional config file with post processing info
    if args.save_config:
        out_config = os.path.join(config[KEY_OUTDIR], OUTPUT_CONFIG)
        with open(out_config, 'w') as f:
            yaml.dump(config, f)

    # runs the preprocessing snakemake
    status = misc.run_snakemake(config,preprocessing_snakefile,config)

    if status: # translate "success" into shell exit code of 0
        with open(os.path.join(config[KEY_TEMPDIR],PREPROCESSING_CONFIG),"r") as f:
            preprocessing_config = yaml.safe_load(f)
        
        status = misc.run_snakemake(preprocessing_config,snakefile,config)

        if status: 

            config[KEY_SAMPLE_SEQS]=os.path.join(config[KEY_OUTDIR],"published_data",SAMPLE_SEQS)

            # initiate phylo 
            if config[KEY_RUN_PHYLO]:

                phylo_outdir = os.path.join(config[KEY_OUTDIR],"phylogenetics")
                if not os.path.exists(phylo_outdir):
                    os.mkdir(phylo_outdir)
                
                # figures out how many trees need building, and what sequences are available to go into them
                # creates the annotations files for the phylogenies too
                seq_clusters,tree_annotations = phylo_functions.get_seqs_and_clusters(config[KEY_SAMPLE_SEQS],
                                                                    config[KEY_SUPPLEMENTARY_SEQUENCES],
                                                                    config[KEY_REFERENCE_SEQUENCES],
                                                                    config[KEY_OUTGROUP_SEQUENCES],
                                                                    config[KEY_BARCODES_CSV],
                                                                    config[KEY_SUPPLEMENTARY_METADATA],
                                                                    phylo_outdir,
                                                                    config)
                config[KEY_CLUSTERS] = seq_clusters
                config[KEY_ANNOTATIONS] = os.path.join(config[KEY_OUTDIR],"phylogenetics","annotations.csv")
                config[KEY_TREE_ANNOTATIONS] = tree_annotations

                #run phylo snakemake
                print(green("Initializing phylo pipeline."))
                status = misc.run_snakemake(config,phylo_snakefile,config)
                
            # get the inputs for making the overall report
            report =os.path.join(config[KEY_OUTDIR],OUTPUT_REPORT)
            summary_csv=os.path.join(config[KEY_TEMPDIR],PREPROCESSING_SUMMARY)
            composition_csv=os.path.join(config[KEY_TEMPDIR],SAMPLE_COMPOSITION)
            detailed_csv = os.path.join(config[KEY_OUTDIR],"detailed_run_report.csv")

            # amalgamate all info to make final html report
            make_output_report(report,
                                config[KEY_BARCODES_CSV],
                                summary_csv,
                                composition_csv,
                                config[KEY_SAMPLE_SEQS],
                                detailed_csv,
                                config[KEY_ANNOTATIONS],
                                config)

            if config[KEY_UPDATE_LOCAL_DATABASE]:
                new_db_seqs = os.path.join(config[KEY_SUPPLEMENTARY_DATADIR],f"{config[KEY_RUNNAME]}.fasta")
                new_db_metadata = os.path.join(config[KEY_SUPPLEMENTARY_DATADIR],f"{config[KEY_RUNNAME]}.csv")
                phylo_functions.update_local_database(config[KEY_SAMPLE_SEQS],detailed_csv,new_db_seqs,new_db_metadata,config)

            for r,d,f in os.walk(os.path.join(config[KEY_OUTDIR],"published_data")):
                for fn in f:
                    if not os.path.getsize(os.path.join(r,fn)):
                        os.remove(os.path.join(r,fn))

            if config[KEY_ARCHIVE_FASTQ]:

                shutil.copytree(config[KEY_READDIR],config[KEY_ARCHIVEDIR],dirs_exist_ok=True)
                print(green("Archiving fastq_pass to ") + f"{config[KEY_ARCHIVEDIR]}.")

            return 0
        return 1
    return 1

if __name__ == '__main__':
    main()
