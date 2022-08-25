# I GROUP KEYS
KEY_BARCODES_CSV = "barcodes_csv"
KEY_READDIR = "readdir"
KEY_BARCODES = "barcodes"
KEY_SAMPLES = "samples"
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
KEY_MEDAKA_MODEL = "medaka_model"

KEY_BARCODE = "barcode"
KEY_SAMPLE = "sample"
KEY_DATE="date"
KEY_REFERENCE = "reference"
KEY_NUM_READS = "num_reads"
KEY_PERCENT = "percent_of_sample"
KEY_REFERENCE_GROUP = "reference_group"
KEY_VIRUS_PRESENT = "virus_present"
KEY_TAXA_OUTDIR = "taxa_outdir"
KEY_HITS = "hits"

# REPORT KEYS
KEY_SNIPIT_SVG="snipit_svg"
KEY_SNP_SITES = "snp_sites"
KEY_SITE = "site"
KEY_INDEL_SITES = "indel_sites"
KEY_MASKED_SITES = "masked_sites"
KEY_SUMMARY_DATA="summary_data"
KEY_VARIATION_INFO="variation_info"
KEY_COOCCURRENCE_INFO="cooccurrence_info"
KEY_POSITIVE="positive"
VALUE_POSITIVE="positive"
KEY_NEGATIVE = "negative"
VALUE_NEGATIVE = "negative"
KEY_SUMMARY_TABLE="summary_table"
KEY_COMPOSITION_TABLE="composition_table"
KEY_COMPOSITION_TABLE_HEADER="composition_table_header"
KEY_SUMMARY_TABLE_HEADER="summary_table_header"
KEY_CONTROL_STATUS="control_status"

# MISC KEYS
KEY_USERNAME="username"
KEY_INSTITUTE="institute"
KEY_RUN_NAME="run_name"
KEY_VERBOSE="verbose"
KEY_THREADS="threads"
KEY_LOG_API="log_api"
KEY_LOG_STRING="log_string"
KEY_QUIET="quiet"
KEY_COLOUR_MAP = "colour_map"
KEY_COLOUR_THEME="colour_theme"
KEY_FASTA="fasta"

KEY_REPORT_TEMPLATE = "report_template"
KEY_BARCODE_REPORT_TEMPLATE = "barcode_report_template"

KEY_LANGUAGE = "language"
VALUE_LANGUAGE = "English"

RESOURCE_KEY_FILENAME="filename"
RESOURCE_KEY_DIRECTORY="directory"
RESOURCE_KEY="KEY"

KEY_SUMMARY_HEADERS = "report_summary_headers"

# default values for config dict
VALUE_OUTPUT_PREFIX = "analysis"
VALUE_SUMMARY_HEADERS = ["taxon","sites","haplotype","num_reads","make_cns"]
VALUE_REFERENCES_FOR_CNS = ["Sabin1-related","Sabin2-related","Sabin3-related","WPV1"]

VALUE_ANALYSIS_MODE = "vp1"
VALUE_ANALYSIS_MODE_VP1 = "vp1"
VALUE_ANALYSIS_MODE_WG_2TILE = "wg_2tile"

READ_LENGTH_DEFAULT_VP1 = [1000,1300]
READ_LENGTH_DEFAULT_WG_2TILE = [3400,5200]

VALUE_DEFAULT_MEDAKA_MODEL="r941_min_high_g360"

VALUE_MIN_READS = 50
VALUE_MIN_PCENT = 10

# vdpv call thresholds

CALL_THRESHOLD_DICT = {
    "Sabin1-related":10,
    "Sabin2-related":8,
    "Sabin3-related":10
}


# report defaults

VALUE_RUN_NAME="Nanopore sequencing"
VALUE_COLOUR_MAP=["#e68781","#476970","#f5eece"]
VALUE_COLOUR_THEME="#e68781"

#file headers
VARIANT_CALLS_HEADER_FIELDS = ["barcode","reference","variant_count","variants"]
SAMPLE_SUMMARY_TABLE_HEADER_FIELDS = ["sample","barcode","Sample call","reference_group","Number of mutations"]
SAMPLE_HIT_HEADER_FIELDS = ["barcode","reference","reference_group","num_reads","percent_of_sample"]

SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_VP1 = ["sample","barcode","Sabin1-related","Sabin2-related","Sabin3-related",
                                "WPV1","WPV2","WPV3","NonPolioEV","unmapped"]
SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_WG = ["sample","barcode","Sabin1-related","Sabin2-related","Sabin3-related","nOPV2",
                                "WPV1","WPV2","WPV3","NonPolioEV","unmapped"]


# file names
OUTPUT_REPORT = "report.html"
SAMPLE_COMPOSITION = "sample_composition.csv"
PREPROCESSING_SUMMARY = "preprocessing_summary.csv"
PREPROCESSING_CONFIG = "preprocessing_config.yaml"
SAMPLE_SEQS = "vp1_sequences.fasta"
REFERENCE_SEQUENCES_FILE_WG = "references.wg.fasta"
REFERENCE_SEQUENCES_FILE_VP1 = "references.vp1.fasta"

# DEPENDENCIES AND RESOURCES TO CHECK
valid_analysis_modes = ["vp1","wg_2tile"]
dependency_list = ["minimap2","snakemake","medaka","racon"]
module_list = ["mako","Bio"]

ENGLISH_RESOURCES = [{RESOURCE_KEY:"report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"english_report.mako"},
        {RESOURCE_KEY:"barcode_report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"english_barcode_report.mako"}
        ]
FRENCH_RESOURCES = [
        {RESOURCE_KEY:"report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"french_report.mako"},
        {RESOURCE_KEY:"barcode_report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"french_barcode_report.mako"}
    ]
