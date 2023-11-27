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
KEY_SAMPLE_TYPE = "sample_type"
KEY_ANALYSIS_MODE = "analysis_mode"
KEY_REFERENCES_FOR_CNS = "references_for_cns"
KEY_MIN_READ_LENGTH = "min_read_length"
KEY_MAX_READ_LENGTH = "max_read_length"
KEY_MIN_READS = "min_read_depth"
KEY_MIN_PCENT = "min_read_pcent"
KEY_MIN_MAP_QUALITY = "min_map_quality"
KEY_MIN_ALN_BLOCK = "min_aln_block"
KEY_REFERENCE_SEQUENCES = "reference_sequences"
KEY_REFERENCE_MATCH_FIELD = "display_name"
KEY_MEDAKA_MODEL = "medaka_model"
KEY_PRIMER_LENGTH = "primer_length"

KEY_BARCODE = "barcode"
KEY_SAMPLE = "sample"
KEY_NAME = "name"
KEY_DATE="sample_date"
KEY_EPID="EPID"
KEY_VARIANTS = "variants"
KEY_VARIANT_COUNT = "variant_count"
KEY_CALL = "call"
KEY_WELL = "well"
KEY_ALL_METADATA="all_metadata_to_header"
KEY_REFERENCE = "reference"
KEY_NUM_READS = "num_reads"
KEY_PERCENT = "percent_of_sample"
KEY_REFERENCE_GROUP = "reference_group"
KEY_SOURCE = "source"
KEY_VIRUS_PRESENT = "virus_present"
KEY_TAXA_OUTDIR = "taxa_outdir"
KEY_HITS = "hits"

KEY_CNS_STEM = "CNS_ID"
VALUE_CNS_STEM = "SEQ"

# PHYLO KEYS

KEY_RUN_PHYLO = "run_phylo"
KEY_UPDATE_LOCAL_DATABASE = "update_local_database"
KEY_SUPPLEMENTARY_DATADIR = "supplementary_datadir"
KEY_SUPPLEMENTARY_SEQUENCES = "supplementary_sequences"
KEY_SUPPLEMENTARY_METADATA = "supplementary_metadata"
KEY_SUPPLEMENTARY_METADATA_ID_COLUMN = "supplementary_metadata_id_column"
KEY_CLUSTERS= "clusters"
KEY_OUTGROUP_SEQUENCES = "outgroup_sequences"
KEY_ANNOTATIONS = "annotations"
KEY_TREE_ANNOTATIONS = "tree_annotations"
KEY_LOCATION = "location"
KEY_LOCAL_DATABASE_THRESHOLD = "local_database_threshold"

KEY_SAMPLE_SEQS = "sample_seqs"
KEY_PHYLO_METADATA_COLUMNS = "phylo_metadata_columns"
KEY_SUPPLEMENTARY_METADATA_COLUMNS = "supplementary_metadata_columns"

VALUE_PHYLO_METADATA_COLUMNS = [KEY_CALL,KEY_DATE,KEY_EPID]
VALUE_SUPPLEMENTARY_METADATA_COLUMNS = [KEY_LOCATION,"lineage"]
VALUE_SUPPLEMENTARY_METADATA_ID_COLUMN = "sequence_name"
VALUE_LOCAL_DATABASE_THRESHOLD = 6

# HAPLO KEYS
KEY_HAPLOTYPE_SAMPLE_SIZE = "haplotype_sample_size"

# REPORT KEYS
KEY_SNIPIT_SVG="snipit_svg"
KEY_SNP_SITES = "snp_sites"
KEY_SITE = "site"
KEY_INDEL_SITES = "indel_sites"
KEY_MASKED_SITES = "masked_sites"
KEY_SUMMARY_DATA="summary_data"
KEY_VARIATION_INFO="variation_info"
KEY_COOCCURRENCE_INFO="cooccurrence_info"
KEY_POSITIVE="positive_control"
KEY_NEGATIVE = "negative_control"
KEY_SUMMARY_TABLE="summary_table"
KEY_COMPOSITION_TABLE="composition_table"
KEY_COMPOSITION_TABLE_HEADER="composition_table_header"
KEY_SUMMARY_TABLE_HEADER="summary_table_header"
KEY_DETAILED_TABLE_HEADER = "detailed_table_header"
KEY_PHYLO_HEADER = "phylo_header"
KEY_CONTROL_STATUS="control_status"
KEY_ORIENTATION="orientation"
KEY_CONFIGURATION_TABLE_FIELDS = "configuration_table_fields"

# MISC KEYS
KEY_USERNAME="username"
KEY_INSTITUTE="institute"
KEY_RUNNAME="runname"
KEY_NOTES="notes"
KEY_VERBOSE="verbose"
KEY_THREADS="threads"
KEY_LOG_API="log_api"
KEY_LOG_STRING="log_string"
KEY_QUIET="quiet"
KEY_COLOUR_MAP = "colour_map"
KEY_COLOUR_THEME="colour_theme"
KEY_FASTA="fasta"
KEY_TODAY = "today"

KEY_REPORT_TEMPLATE = "report_template"
KEY_BARCODE_REPORT_TEMPLATE = "barcode_report_template"

KEY_LANGUAGE = "language"

RESOURCE_KEY_FILENAME="filename"
RESOURCE_KEY_DIRECTORY="directory"
RESOURCE_KEY="KEY"

KEY_SUMMARY_HEADERS = "report_summary_headers"

# default values for config dict

VALUE_LANGUAGE = "English"
VALUE_POSITIVE="positive"
VALUE_NEGATIVE = "negative"

VALUE_OUTPUT_PREFIX = "analysis"
VALUE_SUMMARY_HEADERS = ["taxon","sites","haplotype","num_reads","make_cns"]
VALUE_REFERENCES_FOR_CNS = ["Sabin1-related","Sabin2-related","Sabin3-related","WPV1"]
VALUE_PHYLO_HEADER = ["name","sample","barcode","source","reference_group","call"]

VALUE_SAMPLE_TYPE = "stool"
VALUE_ANALYSIS_MODE = "vp1"
VALUE_ANALYSIS_MODE_VP1 = "vp1"
VALUE_ANALYSIS_MODE_WG = "wg"
VALUE_ANALYSIS_MODE_PANEV = "panev"

READ_LENGTH_DICT = {
    VALUE_ANALYSIS_MODE_WG:[3400,5200],
    VALUE_ANALYSIS_MODE_VP1: [1000,1300],
    VALUE_ANALYSIS_MODE_PANEV:[3000,4500]
}

# READ_LENGTH_DEFAULT_VP1 = [1000,1300]
# READ_LENGTH_DEFAULT_PANEV = [3000,4500]
# READ_LENGTH_DEFAULT_WG = [3400,5200]
VALUE_PRIMER_LENGTH = 30
VALUE_MIN_MAP_QUALITY = 50
VALUE_DEFAULT_MEDAKA_MODEL="r941_min_hac_variant_g507"

VALUE_MIN_READS = 50
VALUE_MIN_PCENT = 10

# vdpv call thresholds

CALL_THRESHOLD_DICT = {
    "Sabin1-related":10,
    "Sabin2-related":6,
    "Sabin3-related":10
}


# ref group default
KEY_DISPLAY_NAME = "display_name"

# report defaults
VALUE_ORIENTATION="vertical"
VALID_ORIENTATION=["vertical","horizontal"]
VALUE_RUNNAME="polioDDNS"
VALUE_COLOUR_MAP=["#e68781","#476970","#fbedac"]
VALUE_COLOUR_THEME="#e68781"

# phylo defaults
VALUE_TREE_ANNOTATIONS="source "

#file headers
VARIANT_CALLS_HEADER_FIELDS = ["barcode","reference","variant_count","variants"]
SAMPLE_SUMMARY_TABLE_HEADER_FIELDS = ["sample","barcode","Sample classification","reference_group","Number of mutations"]
SAMPLE_HIT_HEADER_FIELDS = ["barcode","reference","reference_group","num_reads","percent_of_sample"]

SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_VP1 = ["sample","barcode","Sabin1-related","Sabin2-related","Sabin3-related",
                                "WPV1","WPV2","WPV3","NonPolioEV","unmapped"]


SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_WG = ["sample","barcode","Sabin1-related","Sabin2-related","Sabin3-related","nOPV2",
                                "WPV1","WPV2","WPV3","NonPolioEV","unmapped"]

DETAILED_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_VP1 = [
                    "Sabin1-related|closest_reference","Sabin1-related|num_reads","Sabin1-related|nt_diff_from_reference","Sabin1-related|pcent_match","Sabin1-related|classification",
                    "Sabin2-related|closest_reference","Sabin2-related|num_reads","Sabin2-related|nt_diff_from_reference","Sabin2-related|pcent_match","Sabin2-related|classification",
                    "Sabin3-related|closest_reference","Sabin3-related|num_reads","Sabin3-related|nt_diff_from_reference","Sabin3-related|pcent_match","Sabin3-related|classification",
                    "WPV1|closest_reference","WPV1|num_reads","WPV1|nt_diff_from_reference","WPV1|pcent_match","WPV1|classification",
                    "WPV2|closest_reference","WPV2|num_reads","WPV2|nt_diff_from_reference","WPV2|pcent_match","WPV2|classification",
                    "WPV3|closest_reference","WPV3|num_reads","WPV3|nt_diff_from_reference","WPV3|pcent_match","WPV3|classification",
                    "NonPolioEV|closest_reference","NonPolioEV|num_reads","NonPolioEV|nt_diff_from_reference","NonPolioEV|pcent_match","NonPolioEV|classification","comments"]

DETAILED_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS_WG = [
                    "Sabin1-related|closest_reference","Sabin1-related|num_reads","Sabin1-related|nt_diff_from_reference","Sabin1-related|pcent_match","Sabin1-related|classification",
                    "Sabin2-related|closest_reference","Sabin2-related|num_reads","Sabin2-related|nt_diff_from_reference","Sabin2-related|pcent_match","Sabin2-related|classification",
                    "Sabin3-related|closest_reference","Sabin3-related|num_reads","Sabin3-related|nt_diff_from_reference","Sabin3-related|pcent_match","Sabin3-related|classification",
                    "nOPV2|closest_reference","nOPV2|num_reads","nOPV2|nt_diff_from_reference","nOPV2|pcent_match","nOPV2|classification",
                    "WPV1|closest_reference","WPV1|num_reads","WPV1|nt_diff_from_reference","WPV1|pcent_match","WPV1|classification",
                    "WPV2|closest_reference","WPV2|num_reads","WPV2|nt_diff_from_reference","WPV2|pcent_match","WPV2|classification",
                    "WPV3|closest_reference","WPV3|num_reads","WPV3|nt_diff_from_reference","WPV3|pcent_match","WPV3|classification",
                    "NonPolioEV|closest_reference","NonPolioEV|num_reads","NonPolioEV|nt_diff_from_reference","NonPolioEV|pcent_match","NonPolioEV|classification","comments"]

VALUE_CONFIGURATION_TABLE_FIELDS = [
                    KEY_MIN_READ_LENGTH,KEY_MAX_READ_LENGTH,KEY_MIN_MAP_QUALITY,KEY_MIN_READS,KEY_MIN_ALN_BLOCK,
                    KEY_MIN_PCENT,KEY_MEDAKA_MODEL,KEY_PRIMER_LENGTH,KEY_ANALYSIS_MODE,KEY_RUN_PHYLO,
                    KEY_READDIR
                    ]

# file names
OUTPUT_REPORT = "report.html"
SAMPLE_COMPOSITION = "sample_composition.csv"
PREPROCESSING_SUMMARY = "preprocessing_summary.csv"
PREPROCESSING_CONFIG = "preprocessing_config.yaml"
PHYLO_CONFIG = "phylo_config.yaml"
OUTPUT_CONFIG = "all_config.yaml"
SAMPLE_SEQS = "vp1_sequences.fasta"
REFERENCE_SEQUENCES_FILE_WG = "references.wg.fasta"
REFERENCE_SEQUENCES_FILE_VP1 = "references.vp1.fasta"
OUTGROUP_SEQUENCES_FILE_WG = "outgroups.wg.fasta"
OUTGROUP_SEQUENCES_FILE_VP1 = "outgroups.vp1.fasta"

# DEPENDENCIES AND RESOURCES TO CHECK
valid_analysis_modes = ["vp1","panev","wg"]
valid_sample_types = ["stool","environmental"]

DEPENDENCY_LIST = ["minimap2","snakemake","medaka",]
MODULE_LIST = ["mako","Bio"]
PHYLO_DEPENDENCY_LIST = ["iqtree","mafft","jclusterfunk"]
PHYLO_MODULE_LIST = []

ENGLISH_RESOURCES = [{RESOURCE_KEY:"report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"report.mako"},
        {RESOURCE_KEY:"barcode_report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"barcode_report.mako"}
        ]
FRENCH_RESOURCES = [
        {RESOURCE_KEY:"report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"report.mako"},
        {RESOURCE_KEY:"barcode_report_template",
        RESOURCE_KEY_DIRECTORY:"data",
        RESOURCE_KEY_FILENAME:"barcode_report.mako"}
    ]
