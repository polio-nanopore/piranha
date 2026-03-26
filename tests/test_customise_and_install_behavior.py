import csv
import sys
from pathlib import Path

import pytest

from piranha.input_parsing import customise_run
from piranha.utils import data_install_checks
from piranha.utils.config import (
    KEY_ANALYSIS_MODE,
    KEY_LANGUAGE,
    KEY_OUTGROUP_SEQUENCES,
    KEY_REFERENCE_SEQUENCES,
)


def test_look_for_basecalled_reads_uses_cli_path(tmp_path: Path, monkeypatch):
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()

    monkeypatch.setattr(customise_run, "sys", sys, raising=False)

    config = {}
    customise_run.look_for_basecalled_reads("reads", str(tmp_path), config)
    assert config["read_path"] == str(reads_dir)


def test_look_for_basecalled_reads_uses_config_path(tmp_path: Path, monkeypatch):
    reads_dir = tmp_path / "reads"
    reads_dir.mkdir()
    (reads_dir / "a.fastq").write_text("@r1\nATGC\n+\n!!!!\n")

    monkeypatch.setattr(customise_run, "sys", sys, raising=False)

    config = {"read_path": "reads", "path_to_config": str(tmp_path)}
    customise_run.look_for_basecalled_reads(None, str(tmp_path), config)
    assert config["read_path"] == str(reads_dir)


def test_look_for_barcodes_csv_parses_values(tmp_path: Path, monkeypatch):
    barcode_csv = tmp_path / "barcodes.csv"
    with barcode_csv.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=["barcode"])
        writer.writeheader()
        writer.writerow({"barcode": "NB01"})
        writer.writerow({"barcode": "BC02"})

    monkeypatch.setattr(customise_run, "sys", sys, raising=False)
    monkeypatch.setattr(customise_run, "csv", csv, raising=False)

    config = {}
    customise_run.look_for_barcodes_csv("barcodes.csv", str(tmp_path), config)
    assert config["barcodes"] == "NB01,BC02"


def test_look_for_barcodes_csv_with_no_input_sets_empty(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(customise_run, "sys", sys, raising=False)
    monkeypatch.setattr(customise_run, "csv", csv, raising=False)

    config = {}
    customise_run.look_for_barcodes_csv(None, str(tmp_path), config)
    assert config["barcodes_csv"] == ""
    assert config["barcodes"] == ""


def test_package_data_check_sets_existing_resource_path():
    config = {}
    data_install_checks.package_data_check("references.vp1.fasta", "data", KEY_REFERENCE_SEQUENCES, config)
    assert config[KEY_REFERENCE_SEQUENCES].endswith("references.vp1.fasta")


def test_get_snakefile_finds_existing_file():
    package_dir = Path(__file__).resolve().parents[1] / "piranha"
    snakefile = data_install_checks.get_snakefile(str(package_dir), "main")
    assert snakefile.endswith("piranha_main.smk")


def test_get_references_and_outgroups_choose_by_mode(monkeypatch):
    calls = []

    def fake_package_data_check(filename, directory, key, config):
        calls.append((filename, key))
        config[key] = f"/{directory}/{filename}"

    monkeypatch.setattr(data_install_checks, "package_data_check", fake_package_data_check)

    config = {KEY_ANALYSIS_MODE: "vp1"}
    data_install_checks.get_references(config)
    data_install_checks.get_outgroups(config)

    assert config[KEY_REFERENCE_SEQUENCES].endswith("references.vp1.fasta")
    assert config[KEY_OUTGROUP_SEQUENCES].endswith("outgroups.vp1.fasta")
    assert calls


def test_check_install_rejects_invalid_language():
    config = {KEY_LANGUAGE: "English", KEY_ANALYSIS_MODE: "vp1"}
    with pytest.raises(SystemExit):
        data_install_checks.check_install("Klingon", config)
