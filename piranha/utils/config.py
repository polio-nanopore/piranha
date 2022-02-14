# I GROUP KEYS
KEY_BARCODES_CSV = "barcodes_csv"
KEY_READDIR = "readdir"
KEY_BARCODES = "barcodes"
KEY_BARCODES_TO_SAMPLE = "barcodes_to_samples"
KEY_RUNID = "runid"

# OUTPUT AND DIRECTORIES
KEY_INPUT_PATH = "input_path"
KEY_CWD="cwd"
KEY_OUTDIR = 'outdir'
KEY_PUBLISHDIR = "publishdir"
KEY_TEMPDIR = "tempdir"
KEY_NO_TEMP = "no_temp"
KEY_OUTPUT_PREFIX="output_prefix"
KEY_DATESTAMP="datestamp"
KEY_OVERWRITE="overwrite"

KEY_OUTPUT_REPORT="output_report"

# ANALYSIS OPTION KEYS
KEY_ANALYSIS_MODE = "analysis_mode"
KEY_REFERENCES_FOR_CNS = "references_for_cns"
KEY_MIN_READ_LENGTH = "min_read_length"
KEY_MAX_READ_LENGTH = "max_read_length"
KEY_MIN_READS = "min_read_depth"
KEY_MIN_PCENT = "min_read_pcent"
KEY_REFERENCE_SEQUENCES = "reference_sequences"

KEY_BARCODE = "barcode"
KEY_TAXA_PRESENT = "taxa_present"
KEY_TAXA_OUTDIR = "taxa_outdir"
KEY_HITS = "hits"

# MISC KEYS
KEY_VERBOSE="verbose"
KEY_THREADS="threads"
KEY_LOG_API="log_api"
KEY_LOG_STRING="log_string"
KEY_QUIET="quiet"
KEY_COLOUR_MAP = "colour_map"
KEY_COLOUR_THEME="colour_theme"

KEY_REPORT_TEMPLATE = "report_template"

RESOURCE_KEY_FILENAME="filename"
RESOURCE_KEY_DIRECTORY="directory"
RESOURCE_KEY="KEY"

KEY_SUMMARY_HEADERS = "report_summary_headers"


# DEPENDENCIES AND RESOURCES TO CHECK
valid_analysis_modes = ["stool","environmental"]
dependency_list = ["minimap2","snakemake","medaka","racon"]
module_list = ["mako","Bio"]

resources = [
        {RESOURCE_KEY:"reference_sequences",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"references.VP1.fasta"},
        {RESOURCE_KEY:"report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"report_template.mako"}
    ]