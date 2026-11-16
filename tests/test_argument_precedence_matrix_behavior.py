import ast
from pathlib import Path

import pytest

from piranha.input_parsing import analysis_arg_parsing, initialising
from piranha.utils import misc
from piranha.utils.config import (
	KEY_ALL_METADATA,
	KEY_ANALYSIS_MODE,
	KEY_ARCHIVE_FASTQ,
	KEY_ARCHIVEDIR,
	KEY_BARCODES_CSV,
	KEY_CWD,
	KEY_DATESTAMP,
	KEY_EPI_CSV,
	KEY_HAPLOTYPE_SAMPLE_SIZE,
	KEY_INSTITUTE,
	KEY_LANGUAGE,
	KEY_LOCAL_DATABASE_THRESHOLD,
	KEY_MAX_HAPLOTYPES,
	KEY_MAX_READ_LENGTH,
	KEY_MIN_ALLELE_FREQUENCY,
	KEY_MIN_ALN_BLOCK,
	KEY_MIN_BASE_QUALITY,
	KEY_MIN_HAPLOTYPE_DEPTH,
	KEY_MIN_HAPLOTYPE_DISTANCE,
	KEY_MIN_HIGH_QUALITY_BASE,
	KEY_MIN_MAP_QUALITY,
	KEY_MIN_PCENT,
	KEY_MIN_READ_LENGTH,
	KEY_MIN_READS,
	KEY_MINIMAP2_OPTIONS,
	KEY_NEGATIVE,
	KEY_NO_TEMP,
	KEY_NOTES,
	KEY_ORIENTATION,
	KEY_OUTDIR,
	KEY_OUTPUT_PREFIX,
	KEY_OVERWRITE,
	KEY_PHYLO_METADATA_COLUMNS,
	KEY_POSITIVE,
	KEY_POSITIVE_REFERENCES,
	KEY_PRIMER_LENGTH,
	KEY_PUBLISHDIR,
	KEY_READDIR,
	KEY_REFERENCE_GROUP_FIELD,
	KEY_REFERENCE_SEQUENCES,
	KEY_RUN_HAPLOTYPING,
	KEY_RUN_PHYLO,
	KEY_RUNNAME,
	KEY_SAMPLE_TYPE,
	KEY_SAVE_CONFIG,
	KEY_SUPPLEMENTARY_DATADIR,
	KEY_SUPPLEMENTARY_METADATA_COLUMNS,
	KEY_SUPPLEMENTARY_METADATA_ID_COLUMN,
	KEY_TEMPDIR,
	KEY_THREADS,
	KEY_UPDATE_LOCAL_DATABASE,
	KEY_USERNAME,
	KEY_VERBOSE,
)


ADD_ARG_CASES = [
	("sample_type", KEY_SAMPLE_TYPE, "environmental", "stool"),
	("analysis_mode", KEY_ANALYSIS_MODE, "wg", "vp1"),
	("reference_group_field", KEY_REFERENCE_GROUP_FIELD, "yaml_group", "cli_group"),
	("negative_control", KEY_NEGATIVE, "yaml_negative", "cli_negative"),
	("positive_control", KEY_POSITIVE, "yaml_positive", "cli_positive"),
	("positive_references", KEY_POSITIVE_REFERENCES, "RefY", "RefC"),
	("min_map_quality", KEY_MIN_MAP_QUALITY, 7, 0),
	("min_high_quality_base", KEY_MIN_HIGH_QUALITY_BASE, 12.0, 0.0),
	("min_base_quality", KEY_MIN_BASE_QUALITY, 11, 0),
	("min_read_length", KEY_MIN_READ_LENGTH, 1100, 1200),
	("max_read_length", KEY_MAX_READ_LENGTH, 1400, 1500),
	("min_read_depth", KEY_MIN_READS, 99, 0),
	("min_read_pcent", KEY_MIN_PCENT, 12.5, 0.0),
	("min_aln_block", KEY_MIN_ALN_BLOCK, 0.7, 0.8),
	("primer_length", KEY_PRIMER_LENGTH, 40, 20),
	("run_haplotyping", KEY_RUN_HAPLOTYPING, False, True),
	("haplotype_sample_size", KEY_HAPLOTYPE_SAMPLE_SIZE, 2000, 3000),
	("min_allele_frequency", KEY_MIN_ALLELE_FREQUENCY, 0.2, 0.1),
	("max_haplotypes", KEY_MAX_HAPLOTYPES, 5, 6),
	("min_haplotype_distance", KEY_MIN_HAPLOTYPE_DISTANCE, 3, 4),
	("min_haplotype_depth", KEY_MIN_HAPLOTYPE_DEPTH, 21, 30),
	("run_phylo", KEY_RUN_PHYLO, False, True),
	("update_local_database", KEY_UPDATE_LOCAL_DATABASE, False, True),
	("local_database_threshold", KEY_LOCAL_DATABASE_THRESHOLD, 9.0, 6.0),
	("output_prefix", KEY_OUTPUT_PREFIX, "yaml_output", "cli_output"),
	("datestamp", KEY_DATESTAMP, False, True),
	("overwrite", KEY_OVERWRITE, False, True),
	("no_temp", KEY_NO_TEMP, False, True),
	("all_metadata_to_header", KEY_ALL_METADATA, False, True),
	("language", KEY_LANGUAGE, "French", "English"),
	("save_config", KEY_SAVE_CONFIG, False, True),
	("archive_fastq", KEY_ARCHIVE_FASTQ, False, True),
	("runname", KEY_RUNNAME, "yaml_run", "cli_run"),
	("username", KEY_USERNAME, "yaml_user", "cli_user"),
	("institute", KEY_INSTITUTE, "yaml_institute", "cli_institute"),
	("notes", KEY_NOTES, "yaml notes", "cli notes"),
	("orientation", KEY_ORIENTATION, "horizontal", "vertical"),
	("threads", KEY_THREADS, 8, 1),
	("verbose", KEY_VERBOSE, False, True),
	("phylo_metadata_columns", KEY_PHYLO_METADATA_COLUMNS, "yaml_col", "cli_col"),
	("supplementary_metadata_columns", KEY_SUPPLEMENTARY_METADATA_COLUMNS, "yaml_supp", "cli_supp"),
	("supplementary_metadata_id_column", KEY_SUPPLEMENTARY_METADATA_ID_COLUMN, "yaml_id", "cli_id"),
	("readdir", KEY_READDIR, "yaml_reads", "cli_reads"),
]

PATH_CASES = [
	("outdir", KEY_OUTDIR, "yaml_out", "cli_out"),
	("publishdir", KEY_PUBLISHDIR, "yaml_publish", "cli_publish"),
	("tempdir", KEY_TEMPDIR, "yaml_temp", "cli_temp"),
	("archivedir", KEY_ARCHIVEDIR, "yaml_archive", "cli_archive"),
	("supplementary_datadir", KEY_SUPPLEMENTARY_DATADIR, "yaml_suppdir", "cli_suppdir"),
]

FILE_CASES = [
	("barcodes_csv", KEY_BARCODES_CSV, "yaml_barcodes.csv", "cli_barcodes.csv"),
	("epi_csv", KEY_EPI_CSV, "yaml_epi.csv", "cli_epi.csv"),
	("reference_sequences", KEY_REFERENCE_SEQUENCES, "yaml_refs.fasta", "cli_refs.fasta"),
]


def _command_args_in_use() -> set[str]:
	command_file = Path(__file__).resolve().parents[1] / "piranha" / "command.py"
	tree = ast.parse(command_file.read_text())
	names = set()
	for node in ast.walk(tree):
		if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "args":
			names.add(node.attr)
	return names


def test_precedence_matrix_covers_all_command_args():
	covered = {name for name, *_ in ADD_ARG_CASES}
	covered.update({name for name, *_ in PATH_CASES})
	covered.update({name for name, *_ in FILE_CASES})
	covered.update({"minimap2_options"})

	excluded = {"help", "config"}

	used = _command_args_in_use()
	missing = used - covered - excluded
	assert not missing, f"Add precedence tests for command args: {sorted(missing)}"


@pytest.mark.parametrize("_arg_name,key,yaml_value,cli_value", ADD_ARG_CASES)
def test_add_arg_to_config_respects_default_yaml_cli_order(_arg_name, key, yaml_value, cli_value):
	config = initialising.get_defaults()

	config[key] = yaml_value
	misc.add_arg_to_config(key, None, config)
	assert config[key] == yaml_value

	misc.add_arg_to_config(key, cli_value, config)
	assert config[key] == cli_value


@pytest.mark.parametrize("_arg_name,key,yaml_rel,cli_rel", PATH_CASES)
def test_add_path_to_config_respects_default_yaml_cli_order(tmp_path, _arg_name, key, yaml_rel, cli_rel):
	config = {KEY_CWD: str(tmp_path)}

	config[key] = str(tmp_path / yaml_rel)
	misc.add_path_to_config(key, None, config)
	assert config[key] == str(tmp_path / yaml_rel)

	misc.add_path_to_config(key, cli_rel, config)
	assert config[key] == str(tmp_path / cli_rel)


@pytest.mark.parametrize("_arg_name,key,yaml_rel,cli_rel", FILE_CASES)
def test_add_file_to_config_respects_default_yaml_cli_order(tmp_path, _arg_name, key, yaml_rel, cli_rel):
	config = {KEY_CWD: str(tmp_path)}

	config[key] = str(tmp_path / yaml_rel)
	misc.add_file_to_config(key, None, config)
	assert config[key] == str(tmp_path / yaml_rel)

	misc.add_file_to_config(key, cli_rel, config)
	assert config[key] == str(tmp_path / cli_rel)


def test_minimap2_options_precedence_prefers_cli_over_yaml():
	config = {
		KEY_MINIMAP2_OPTIONS: "k=15 w=10",
	}

	parsed = analysis_arg_parsing.minimap2_options_parsing(["k=7", "w=4"], config)
	assert "-k7" in parsed
	assert "-w4" in parsed


def test_minimap2_options_uses_yaml_when_cli_missing():
	config = {
		KEY_MINIMAP2_OPTIONS: "k=15 w=10",
	}

	parsed = analysis_arg_parsing.minimap2_options_parsing(None, config)
	assert "-k15" in parsed
	assert "-w10" in parsed
