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
    KEY_SAVE_CONFIG,
    KEY_TEMPDIR,
    KEY_VERBOSE,
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
    monkeypatch.setattr(command.init, "setup_config_dict", lambda *a, **k: (config, {}))
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


def test_main_honors_yaml_save_config_when_cli_omits_flag(monkeypatch, tmp_path: Path):
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
    config[KEY_SAVE_CONFIG] = True

    (outdir / "published_data").mkdir(exist_ok=True)

    monkeypatch.setattr(command.dependency_checks, "check_dependencies", lambda *a, **k: None)
    monkeypatch.setattr(command.init, "setup_config_dict", lambda *a, **k: (config, {KEY_SAVE_CONFIG: True}))
    monkeypatch.setattr(command.analysis_arg_parsing, "sample_type", lambda *a, **k: None)
    monkeypatch.setattr(command.data_install_checks, "check_install", lambda *a, **k: None)
    monkeypatch.setattr(command.analysis_arg_parsing, "analysis_mode", lambda *a, **k: None)
    monkeypatch.setattr(command.misc, "add_check_valid_arg", lambda *a, **k: None)
    monkeypatch.setattr(command.data_install_checks, "get_snakefile", lambda *a, **k: "Snakefile")
    monkeypatch.setattr(command.analysis_arg_parsing, "analysis_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.analysis_arg_parsing, "minimap2_options_parsing", lambda *a, **k: "-x map-ont")
    monkeypatch.setattr(command.analysis_arg_parsing, "haplo_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "parse_input_group", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "control_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.directory_setup, "output_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.init, "misc_args_to_config", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "phylo_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command, "make_output_report", lambda *a, **k: None)

    captured = {"snake_configfile": None}

    def fake_run_snakemake(snake_configfile, *args, **kwargs):
        captured["snake_configfile"] = snake_configfile
        return False

    monkeypatch.setattr(command.misc, "run_snakemake", fake_run_snakemake)

    exit_code = command.main(["-b", "barcodes.csv", "-i", "reads", "-r", "refs.fasta"])
    assert exit_code == 1
    assert captured["snake_configfile"] == str(outdir / command.OUTPUT_CONFIG)


def test_main_honors_yaml_verbose_when_cli_omits_flag(monkeypatch, tmp_path: Path):
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
    config[KEY_VERBOSE] = True

    (outdir / "published_data").mkdir(exist_ok=True)

    monkeypatch.setattr(command.dependency_checks, "check_dependencies", lambda *a, **k: None)
    monkeypatch.setattr(command.init, "setup_config_dict", lambda *a, **k: (config, {KEY_VERBOSE: True}))
    monkeypatch.setattr(command.analysis_arg_parsing, "sample_type", lambda *a, **k: None)
    monkeypatch.setattr(command.data_install_checks, "check_install", lambda *a, **k: None)
    monkeypatch.setattr(command.analysis_arg_parsing, "analysis_mode", lambda *a, **k: None)
    monkeypatch.setattr(command.misc, "add_check_valid_arg", lambda *a, **k: None)
    monkeypatch.setattr(command.data_install_checks, "get_snakefile", lambda *a, **k: "Snakefile")
    monkeypatch.setattr(command.analysis_arg_parsing, "analysis_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.analysis_arg_parsing, "minimap2_options_parsing", lambda *a, **k: "-x map-ont")
    monkeypatch.setattr(command.analysis_arg_parsing, "haplo_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "parse_input_group", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "control_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.directory_setup, "output_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.init, "misc_args_to_config", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "phylo_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command, "make_output_report", lambda *a, **k: None)

    captured = {"verbose": None}

    def fake_run_snakemake(snake_configfile, snakefile, verbose, cfg, extra_config=None):
        captured["verbose"] = verbose
        return False

    monkeypatch.setattr(command.misc, "run_snakemake", fake_run_snakemake)

    exit_code = command.main(["-b", "barcodes.csv", "-i", "reads", "-r", "refs.fasta"])
    assert exit_code == 1
    assert captured["verbose"] is True


def test_main_end_to_end_precedence_yaml_then_cli(monkeypatch, tmp_path: Path):
    yaml_outdir = tmp_path / "yaml_out"
    yaml_tempdir = tmp_path / "yaml_temp"
    cli_outdir = tmp_path / "cli_out"

    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        "\n".join(
            [
                "analysis_mode: wg",
                "min_read_length: 3650",
                "max_read_length: 5100",
                "threads: 7",
                "verbose: true",
                "save_config: false",
                f"outdir: {yaml_outdir}",
                f"tempdir: {yaml_tempdir}",
            ]
        )
        + "\n"
    )

    monkeypatch.setattr(command.dependency_checks, "check_dependencies", lambda *a, **k: None)
    monkeypatch.setattr(command.data_install_checks, "check_install", lambda *a, **k: None)
    monkeypatch.setattr(command.data_install_checks, "get_snakefile", lambda *a, **k: "Snakefile")
    monkeypatch.setattr(command.input_qc, "parse_input_group", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "control_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command.input_qc, "phylo_group_parsing", lambda *a, **k: None)
    monkeypatch.setattr(command, "make_output_report", lambda *a, **k: None)

    captured = {"snake_configfile": None, "verbose": None, "config": None}

    def fake_run_snakemake(snake_configfile, snakefile, verbose, cfg_dict, extra_config=None):
        captured["snake_configfile"] = snake_configfile
        captured["verbose"] = verbose
        captured["config"] = dict(cfg_dict)
        return False

    monkeypatch.setattr(command.misc, "run_snakemake", fake_run_snakemake)

    exit_code = command.main(
        [
            "-c",
            str(cfg),
            "--threads",
            "2",
            "--save-config",
            "--outdir",
            str(cli_outdir),
        ]
    )

    assert exit_code == 1
    assert captured["config"] is not None

    # CLI should override YAML for explicitly supplied fields.
    assert captured["config"][command.KEY_THREADS] == 2
    assert captured["config"][command.KEY_SAVE_CONFIG] is True
    assert Path(captured["config"][command.KEY_OUTDIR]).parts[-1].startswith("cli_out")

    # YAML should override defaults when CLI omits those fields.
    assert captured["config"][command.KEY_VERBOSE] is True
    assert captured["verbose"] is True
    assert captured["config"][command.KEY_MIN_READ_LENGTH] == 3650
    assert captured["config"][command.KEY_MAX_READ_LENGTH] == 5100
    assert captured["config"][command.KEY_ANALYSIS_MODE] == "wg"

    # save-config true should write OUTPUT_CONFIG to outdir, not tempdir.
    assert captured["snake_configfile"] == str(Path(captured["config"][command.KEY_OUTDIR]) / command.OUTPUT_CONFIG)
