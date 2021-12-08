#!/usr/bin/env python3
from piranha import __version__

from piranha.utils import misc
from piranha.utils import dependency_checks
from piranha.utils import data_install_checks

import piranha.utils.custom_logger as custom_logger
from piranha.utils.log_colours import green,cyan,red
from piranha.utils.config import *

import os
import sys
import argparse
import snakemake

cwd = os.getcwd()
thisdir = os.path.abspath(os.path.dirname(__file__))

def main(sysargs = sys.argv[1:]):
 
    parser = argparse.ArgumentParser(add_help=False,
    description=misc.preamble(__version__),
    usage='''
\tpiranha -c <config.yaml> [options]
\tpiranha -i input.csv [options]
\tpiranha -ids ID1,ID2 [options]
\tpiranha -fm <column=match> [options]\n\n''')

    i_group = parser.add_argument_group('Input options')
    i_group.add_argument('-c',"--config", action="store",help="Input config file in yaml format, all command line arguments can be passed via the config file.", dest="config")
    io_group.add_argument('-i','--read-path',help="Path to the directory containing fastq files",dest="read_path")

    barcode_group = parser.add_argument_group('Barcode options')
    barcode_group.add_argument('-b','--barcodes-csv',help="CSV file describing which barcodes were used on which sample",dest="barcodes_csv")
    barcode_group.add_argument('-k','--barcode-kit',help="Indicates which barcode kit was used. Default: native. Options: native, rapid, pcr, all",dest="barcode_kit")

    demux_group = parser.add_argument_group('Demultiplexing options')
    demux_group.add_argument('--demultiplex',action="store_true",help="Indicates that your reads have not been demultiplexed and will run guppy demultiplex on your provided read directory",dest="demultiplex")
    demux_group.add_argument('--path-to-guppy',action="store",help="Path to guppy_barcoder executable",dest="path_to_guppy")

    run_group = parser.add_argument_group('Analysis options')
    run_group.add_argument("-r","--report",action="store_true",help="Generate markdown report of estimated age")
    
    o_group = parser.add_argument_group('Output options')
    o_group.add_argument('-o','--outdir', action="store",help="Output directory. Default: `piranha-2021-XX-YY`")
    o_group.add_argument('-temp','--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: `$TMPDIR`")
    o_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files. For development/ debugging purposes",dest="no_temp")

    misc_group = parser.add_argument_group('Misc options')
    misc_group.add_argument('-t', '--threads', action='store',dest="threads",type=int,help="Number of threads")
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

    dependency_checks.check_dependencies(dependency_list, module_list)
    if args.mutations:
        dependency_checks.check_scorpio_mutations(args.mutations)
    
    # Initialise config dict
    config = init.setup_config_dict(cwd,args.config)

    # Checks access to package data and grabs the snakefile
    data_install_checks.check_install(config)
    snakefile = data_install_checks.get_snakefile(thisdir)
    
    # Threads and verbosity to config
    init.misc_args_to_config(args.verbose,args.threads,args.piranha_mode,config)
    init.set_up_verbosity(config)

    # Checks access to package data and grabs the snakefile
    data_install_checks.check_install(config)
    
    # Sort out where the query info is coming from, csv or id string, optional fasta seqs.
    # Checks if they're real files, of the right format and that QC args sensible values.

    """
    Configure whether guppy barcoder needs to be run
    """

    qcfunk.look_for_guppy_barcoder(args.demultiplex,args.path_to_guppy,cwd,config)



    snakefile = data_install_checks.get_snakefile(thisdir)

    # sets up the output dir, temp dir, and data output desination
    directory_setup.output_group_parsing(args.outdir, args.output_prefix, args.overwrite, args.datestamp, args.output_data, args.tempdir, args.no_temp, config)
    # ready to run? either verbose snakemake or quiet mode

    if config[KEY_VERBOSE]:
        print(red("\n**** CONFIG ****"))
        for k in sorted(config):
            print(green(f" - {k}: ") + f"{config[k]}")
        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                    workdir=config[KEY_TEMPDIR],config=config, cores=config[KEY_THREADS],lock=False
                                    )
    else:
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=config[KEY_TEMPDIR],
                                    config=config, cores=config[KEY_THREADS],lock=False,quiet=True,log_handler=config[KEY_LOG_API]
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1



if __name__ == '__main__':
    main()
