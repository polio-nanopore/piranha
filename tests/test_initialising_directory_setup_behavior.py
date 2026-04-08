from pathlib import Path
from datetime import date

import pytest

from piranha.input_parsing import directory_setup, initialising
from piranha.utils.config import (
    KEY_ARCHIVE_FASTQ,
    KEY_ARCHIVEDIR,
    KEY_CWD,
    KEY_DATESTAMP,
    KEY_LOG_STRING,
    KEY_NO_TEMP,
    KEY_OUTDIR,
    KEY_OUTPUT_PREFIX,
    KEY_OVERWRITE,
    KEY_QUIET,
    KEY_READDIR,
    KEY_REFERENCE_SEQUENCES,
    KEY_TEMPDIR,
    KEY_VERBOSE,
)


def test_check_configfile_rejects_non_yaml(tmp_path: Path):
    cfg = tmp_path / "config.txt"
    cfg.write_text("a: 1\n")

    with pytest.raises(SystemExit):
        initialising.check_configfile(str(tmp_path), cfg.name)


def test_parse_yaml_file_sets_absolute_paths(tmp_path: Path):
    cfg = tmp_path / "config.yaml"
    cfg.write_text("readdir: reads\nreference_sequences: refs.fasta\n")

    config = initialising.get_defaults()
    initialising.parse_yaml_file(str(cfg), config)

    assert config[KEY_READDIR] == str(tmp_path / "reads")
    assert config[KEY_REFERENCE_SEQUENCES] == str(tmp_path / "refs.fasta")


# def test_set_up_verbosity_sets_quiet_and_log_string():
#     config = {KEY_VERBOSE: False}
#     initialising.set_up_verbosity(config)

#     assert config[KEY_QUIET] is True
#     assert "--quiet --log-handler-script" in config[KEY_LOG_STRING]


def test_datestamped_outdir_increments_when_existing(tmp_path: Path):
    today = date.today().strftime("%Y-%m-%d")
    existing = tmp_path / f"analysis_{today}"
    existing.mkdir()
    config = {
        KEY_CWD: str(tmp_path),
        KEY_OUTPUT_PREFIX: "analysis",
        KEY_DATESTAMP: True,
        KEY_OVERWRITE: False,
        KEY_OUTDIR: str(tmp_path / "analysis"),
    }

    outdir, _ = directory_setup.datestamped_outdir(config)
    assert outdir == f"{existing}_1"


def test_set_up_archivedir_creates_output_subdir(tmp_path: Path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    config = {
        KEY_ARCHIVE_FASTQ: True,
        KEY_OUTDIR: str(outdir),
    }

    directory_setup.set_up_archivedir(config)
    assert (outdir / "fastq_pass").exists()
    assert config[KEY_ARCHIVEDIR] == str(outdir / "fastq_pass")


def test_set_up_tempdir_uses_outdir_when_no_temp(tmp_path: Path):
    outdir = tmp_path / "out"
    outdir.mkdir()
    config = {
        KEY_NO_TEMP: True,
        KEY_OUTDIR: str(outdir),
    }

    directory_setup.set_up_tempdir(config)
    assert config[KEY_TEMPDIR] == str(outdir)


def test_output_group_parsing_creates_outdir_and_tempdir(tmp_path: Path):
    config = initialising.get_defaults()
    config[KEY_CWD] = str(tmp_path)

    directory_setup.output_group_parsing(
        outdir="analysis_out",
        output_prefix="analysis",
        overwrite=False,
        datestamp=False,
        tempdir="",
        no_temp=False,
        archive_fastq=False,
        archivedir="",
        config=config,
    )

    assert Path(config[KEY_OUTDIR]).exists()
    assert Path(config[KEY_TEMPDIR]).exists()
