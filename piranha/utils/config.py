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
KEY_SAMPLE = "sample"
KEY_REFERENCE = "reference"
KEY_NUM_READS = "num_reads"
KEY_PERCENT = "percent_of_sample"
KEY_REFERENCE_GROUP = "reference_group"
KEY_VIRUS_PRESENT = "virus_present"
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

# default values for config dict
VALUE_OUTPUT_PREFIX = "analysis"
VALUE_SUMMARY_HEADERS = ["taxon","sites","haplotype","num_reads","make_cns"]
VALUE_REFERENCES_FOR_CNS = ["Sabin1-related","Sabin2-related","Sabin3-related","WPV1"]
VALUE_ANALYSIS_MODE = "stool"
VALUE_MIN_READ_LENGTH = 1000
VALUE_MAX_READ_LENGTH = 1300
VALUE_MIN_READS = 50
VALUE_MIN_PCENT = 5

VALUE_COLOUR_MAP=["#e68781","#476970","#f5eece"]
VALUE_COLOUR_THEME="#e68781"

#file headers

SAMPLE_HIT_HEADER_FIELDS = ["barcode","reference","reference_group","num_reads","percent_of_sample"]
SAMPLE_SUMMARY_HEADER_FIELDS = ["sample","barcode","Sabin1-related","Sabin2-related","Sabin3-related",
                                "WPV1","WPV2","WPV3","NonPolioEV","unmapped","virus_present"]

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