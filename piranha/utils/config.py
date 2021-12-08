
# OUTPUT AND DIRECTORIES
KEY_INPUT_PATH = "input_path"
KEY_READDIR = "readdir"
KEY_CWD="cwd"
KEY_OUTDIR = 'outdir'
KEY_TEMPDIR = "tempdir"
KEY_NO_TEMP = "no_temp"
KEY_OUTPUT_PREFIX="output_prefix"
KEY_DATESTAMP="datestamp"
KEY_OVERWRITE="overwrite"

# ANALYSIS OPTION KEYS

KEY_MIN_READ_LENGTH = "min_read_length"
KEY_MAX_READ_LENGTH = "max_read_length"
KEY_MIN_READS = "min_read_depth"
KEY_MIN_PCENT = "min_read_pcent"
KEY_REFERENCE_SEQUENCES = "reference_sequences"

# MISC KEYS
KEY_VERBOSE="verbose"
KEY_THREADS="threads"
KEY_LOG_API="log_api"
KEY_LOG_STRING="log_string"
KEY_QUIET="quiet"

# DEPENDENCIES AND RESOURCES TO CHECK
dependency_list = ["gofasta","minimap2","snakemake","iqtree","jclusterfunk"]
module_list = ["mako","Bio"]

resources = [
        {RESOURCE_KEY:"reference_sequence",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"references.fasta"},
        {RESOURCE_KEY:"report_template",
        RESOURCE_KEY_DIRECTORY:"data/report_modules",
        RESOURCE_KEY_FILENAME:"report_template.mako"}
    ]