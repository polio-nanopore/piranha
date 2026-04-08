import pytest
from pathlib import Path
import yaml

from piranha import command
from piranha.input_parsing import initialising
from piranha.utils.config import (
    KEY_ANNOTATIONS,
    KEY_ARCHIVE_FASTQ,
    KEY_BARCODES_CSV,
    KEY_EPI_CSV,
    KEY_OUTDIR,
    KEY_RUN_PHYLO,
    KEY_TEMPDIR,
)


def test_main_without_args_exits_with_help_code_zero():
    with pytest.raises(SystemExit) as exc_info:
        command.main([])

    assert exc_info.value.code == 0


def test_main_help_flag_exits_with_help_code_zero():
    with pytest.raises(SystemExit) as exc_info:
        command.main(["--help"])

    assert exc_info.value.code == 0


def test_main_version_flag_exits_with_code_zero():
    with pytest.raises(SystemExit) as exc_info:
        command.main(["--version"])

    assert exc_info.value.code == 0


def test_main_mocked_success_path_returns_zero(monkeypatch, tmp_path: Path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    tempdir = tmp_path / "temp"
    tempdir.mkdir()

    config = initialising.get_defaults()
    config[KEY_OUTDIR] = str(outdir)
    config[KEY_TEMPDIR] = str(tempdir)
    config[KEY_RUN_PHYLO] = False
    config[KEY_ARCHIVE_FASTQ] = False
    config[KEY_BARCODES_CSV] = ""
    config[KEY_EPI_CSV] = ""
    config[KEY_ANNOTATIONS] = ""

    (outdir / "published_data").mkdir(exist_ok=True)
    (tempdir / command.PREPROCESSING_SUMMARY).write_text("sample,barcode\n")
    (tempdir / command.SAMPLE_COMPOSITION).write_text("sample,barcode\n")

    monkeypatch.setattr(command.dependency_checks, "check_dependencies", lambda *a, **k: None)
    monkeypatch.setattr(command.init, "setup_config_dict", lambda *a, **k: config)
    monkeypatch.setattr(command.analysis_arg_parsing, "sample_type", lambda *a, **k: None)
    monkeypatch.setattr(command.data_install_checks, "check_install", lambda *a, **k: None)
    monkeypatch.setattr(command.analysis_arg_parsing, "analysis_mode", lambda *a, **k: None)
    monkeypatch.setattr(command.misc, "add_check_valid_arg", lambda *a, **k: None)
    monkeypatch.setattr(command.data_install_checks, "get_snakefile", lambda *a, **k: "Snakefile")
    monkeypatch.setattr(command.analysis_arg_parsing, "medaka_options_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.analysis_arg_parsing, "analysis_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.analysis_arg_parsing, "minimap2_options_parsing", lambda *a, **k: "-x map-ont")
    monkeypatch.setattr(command.analysis_arg_parsing, "haplo_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.misc, "add_arg_to_config", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "parse_input_group", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "control_group_parsing", lambda *a, **k: None)

    def fake_output_group_parsing(*args, **kwargs):
        cfg = kwargs["config"] if "config" in kwargs else args[-1]
        cfg[KEY_OUTDIR] = str(outdir)
        cfg[KEY_TEMPDIR] = str(tempdir)

    monkeypatch.setattr(command.directory_setup, "output_group_parsing", fake_output_group_parsing)
    monkeypatch.setattr(command.init, "misc_args_to_config", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "phylo_group_parsing", lambda *a, **k: None)
    # monkeypatch.setattr(command.init, "set_up_verbosity", lambda *a, **k: None)

    calls = {"n": 0}

    def fake_run_snakemake(snake_configfile, snakefile, args_verbose, cfg, extra_config=None):
        calls["n"] += 1
        if calls["n"] == 1:
            with open(Path(cfg[KEY_TEMPDIR]) / command.PREPROCESSING_CONFIG, "w") as handle:
                yaml.safe_dump({"any": "value"}, handle)
            return True
        return True

    monkeypatch.setattr(command.misc, "run_snakemake", fake_run_snakemake)
    monkeypatch.setattr(command.make_output_report, "__call__", command.make_output_report)
    monkeypatch.setattr(command, "make_output_report", lambda *a, **k: None)
    monkeypatch.setattr(command.phylo_functions, "update_local_database", lambda *a, **k: None)

    exit_code = command.main(["-b", "barcodes.csv", "-i", "reads", "-r", "refs.fasta"])
    assert exit_code == 0
