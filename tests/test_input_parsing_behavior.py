import pytest
from pathlib import Path

from piranha.input_parsing import analysis_arg_parsing, input_qc
from piranha.input_parsing import initialising
from piranha.utils.config import (
    KEY_BARCODES,
    KEY_INCLUDE_POSITIVE_REFERENCES,
    KEY_MIN_PCENT,
    KEY_MIN_READS,
    KEY_NEGATIVE,
    KEY_POSITIVE,
    KEY_POSITIVE_REFERENCES,
    KEY_READDIR,
    KEY_REFERENCE_SEQUENCES,
    KEY_RUNID,
    KEY_SAMPLE,
    KEY_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS,
    KEY_DETAILED_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS,
    KEY_SAMPLES,
    KEY_REFERENCE_GROUP_VALUES,
    KEY_RUN_PHYLO,
    KEY_TEMPDIR,
)


def test_define_valid_wells_has_96_wells():
    wells = input_qc.define_valid_wells()
    assert len(wells) == 96
    assert "A01" in wells
    assert "H12" in wells


def test_parse_ref_group_values_extracts_annotation():
    description = "REF001 reference_group=Sabin2-related runname=test"
    assert input_qc.parse_ref_group_values(description, "reference_group") == "Sabin2-related"


def test_pattern_match_for_controls_matches_substring():
    samples = ["neg_ctrl_1", "sampleA", "pos_ctrl_2"]
    assert input_qc.pattern_match_for_controls(samples, "ctrl") == ["neg_ctrl_1", "pos_ctrl_2"]


def test_check_if_int_and_float_convert_values():
    config = {KEY_MIN_READS: "30", KEY_MIN_PCENT: "15.5"}
    analysis_arg_parsing.check_if_int(KEY_MIN_READS, config)
    analysis_arg_parsing.check_if_float(KEY_MIN_PCENT, config)

    assert config[KEY_MIN_READS] == 30
    assert config[KEY_MIN_PCENT] == 15.5


def test_check_if_int_exits_for_invalid_value():
    config = {KEY_MIN_READS: "not-an-int"}
    with pytest.raises(SystemExit):
        analysis_arg_parsing.check_if_int(KEY_MIN_READS, config)


def test_check_if_float_exits_for_invalid_value():
    config = {KEY_MIN_PCENT: "not-a-float"}
    with pytest.raises(SystemExit):
        analysis_arg_parsing.check_if_float(KEY_MIN_PCENT, config)


def test_check_samples_match_passes(tmp_path: Path):
    epi_csv = tmp_path / "epi.csv"
    epi_csv.write_text("sample,field\nS1,x\nS2,y\n")

    input_qc.check_samples_match(str(epi_csv), ["S1", "S2"])


def test_check_samples_match_exits_when_missing_sample(tmp_path: Path):
    epi_csv = tmp_path / "epi.csv"
    epi_csv.write_text("sample,field\nS1,x\n")

    with pytest.raises(SystemExit):
        input_qc.check_samples_match(str(epi_csv), ["S1", "S2"])


def test_parse_read_dir_sets_parent_and_run_id(tmp_path: Path):
    demux = tmp_path / "demultiplexed"
    barcode_dir = demux / "barcode01"
    barcode_dir.mkdir(parents=True)
    (barcode_dir / "reads_runABC_20260326.fastq").write_text("@r1\nATGC\n+\n!!!!\n")

    config = {
        KEY_BARCODES: ["barcode01"],
        KEY_READDIR: "",
    }

    input_qc.parse_read_dir(str(demux), config)
    assert config[KEY_READDIR] == str(demux)
    assert config[KEY_RUNID] == "runABC"


def test_parse_read_dir_exits_when_no_matching_barcodes(tmp_path: Path):
    demux = tmp_path / "demultiplexed"
    barcode_dir = demux / "barcode02"
    barcode_dir.mkdir(parents=True)
    (barcode_dir / "reads_runABC_20260326.fastq").write_text("@r1\nATGC\n+\n!!!!\n")

    config = {
        KEY_BARCODES: ["barcode01"],
        KEY_READDIR: "",
    }

    with pytest.raises(SystemExit):
        input_qc.parse_read_dir(str(demux), config)


def test_parse_barcodes_csv_sets_samples_and_barcodes(tmp_path: Path):
    barcode_csv = tmp_path / "barcodes.csv"
    barcode_csv.write_text("sample,barcode,EPID\nS1,barcode01,E1\nS2,barcode02,E2\n")

    config = {"cwd": str(tmp_path)}
    input_qc.parse_barcodes_csv(str(barcode_csv), config)

    assert config[KEY_BARCODES] == ["barcode01", "barcode02"]
    assert config[KEY_SAMPLES] == ["S1", "S2"]


def test_parse_barcodes_csv_exits_on_duplicate_barcode(tmp_path: Path):
    barcode_csv = tmp_path / "barcodes.csv"
    barcode_csv.write_text("sample,barcode\nS1,barcode01\nS2,barcode01\n")

    config = {"cwd": str(tmp_path)}
    with pytest.raises(SystemExit):
        input_qc.parse_barcodes_csv(str(barcode_csv), config)


def test_control_group_parsing_sets_include_positive_refs(tmp_path: Path):
    refs = tmp_path / "refs.fasta"
    refs.write_text(">RefA reference_group=Sabin2-related\nATGC\n")

    config = {
        KEY_SAMPLES: ["pos_ctrl_1", "neg_ctrl_1", "sampleA"],
        KEY_POSITIVE: "pos_ctrl",
        KEY_NEGATIVE: "neg_ctrl",
        KEY_POSITIVE_REFERENCES: "RefA",
        KEY_INCLUDE_POSITIVE_REFERENCES: False,
        KEY_REFERENCE_SEQUENCES: str(refs),
    }

    input_qc.control_group_parsing("pos_ctrl", "neg_ctrl", "RefA", config)

    assert "pos_ctrl_1" in config[KEY_POSITIVE]
    assert "neg_ctrl_1" in config[KEY_NEGATIVE]
    assert config[KEY_INCLUDE_POSITIVE_REFERENCES] is True


def test_parse_input_group_populates_reference_headers(tmp_path: Path, monkeypatch):
    barcodes = tmp_path / "barcodes.csv"
    refs = tmp_path / "refs.fasta"

    barcodes.write_text("sample,barcode\nS1,barcode01\n")
    refs.write_text(
        ">RefA reference_group=Sabin2-related\nATGC\n"
        ">RefB reference_group=WPV1\nATGC\n"
    )

    config = {
        "cwd": str(tmp_path),
        "path_to_config": str(tmp_path),
        "reference_group_field": "reference_group",
        "epi_csv": "",
    }

    monkeypatch.setattr(input_qc, "parse_read_dir", lambda readdir, cfg: cfg.update({KEY_RUNID: "run1", KEY_READDIR: str(tmp_path / "demux")}))

    input_qc.parse_input_group(
        barcodes_csv=str(barcodes),
        epi_csv="",
        readdir=str(tmp_path / "demux"),
        reference_sequences=str(refs),
        reference_group_field="reference_group",
        config=config,
    )

    assert config[KEY_REFERENCE_GROUP_VALUES] == ["Sabin2-related", "WPV1"]
    assert KEY_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS in config
    assert KEY_DETAILED_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS in config
    assert "Sabin2-related" in config[KEY_SAMPLE_COMPOSITION_TABLE_HEADER_FIELDS]


def test_phylo_group_parsing_true_validates_metadata_columns(tmp_path: Path):
    barcodes = tmp_path / "barcodes.csv"
    barcodes.write_text("sample,barcode\nS1,barcode01\n")
    tempdir = tmp_path / "temp"
    tempdir.mkdir()

    config = initialising.get_defaults()
    config[KEY_TEMPDIR] = str(tempdir)
    config["cwd"] = str(tmp_path)
    config["path_to_config"] = str(tmp_path)

    input_qc.phylo_group_parsing(
        run_phylo_arg=True,
        update_local_database=False,
        supplementary_datadir="",
        phylo_metadata_columns_arg=["sample"],
        barcodes_csv=str(barcodes),
        supplementary_metadata_columns_arg=["country"],
        supplementary_metadata_id_column_arg="name",
        local_database_threshold=10,
        config=config,
    )

    assert config[KEY_RUN_PHYLO] is True
