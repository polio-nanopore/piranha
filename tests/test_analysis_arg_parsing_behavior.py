import pytest

from piranha.input_parsing import analysis_arg_parsing
from piranha.utils.config import (
    ENVIRONMENTAL_DEFAULT_DICT,
    KEY_ANALYSIS_MODE,
    KEY_MAX_READ_LENGTH,
    KEY_MEDAKA_MODEL,
    KEY_MIN_ALLELE_FREQUENCY,
    KEY_MIN_HAPLOTYPE_DEPTH,
    KEY_MIN_HAPLOTYPE_DISTANCE,
    KEY_MIN_READ_LENGTH,
    KEY_RUN_HAPLOTYPING,
    KEY_SAMPLE_TYPE,
    KEY_HAPLOTYPE_SAMPLE_SIZE,
    KEY_MAX_HAPLOTYPES,
)


def test_get_available_medaka_models_parses_stdout(monkeypatch):
    class Result:
        stdout = b"header\nAvailable: r941_min_hac_g507, r1041_e82\n"

    monkeypatch.setattr(analysis_arg_parsing.subprocess, "run", lambda *a, **k: Result())
    models = analysis_arg_parsing.get_available_medaka_models()
    assert models == ["r941_min_hac_g507", "r1041_e82"]


def test_medaka_options_parsing_accepts_valid_model(monkeypatch):
    monkeypatch.setattr(analysis_arg_parsing, "get_available_medaka_models", lambda: ["m1", "m2"])
    config = {KEY_MEDAKA_MODEL: "m1"}

    analysis_arg_parsing.medaka_options_parsing("m2", False, config)
    assert config[KEY_MEDAKA_MODEL] == "m2"


def test_medaka_options_parsing_list_models_exits_zero(monkeypatch):
    monkeypatch.setattr(analysis_arg_parsing, "get_available_medaka_models", lambda: ["m1", "m2"])
    config = {KEY_MEDAKA_MODEL: "m1"}

    with pytest.raises(SystemExit) as exc_info:
        analysis_arg_parsing.medaka_options_parsing(None, True, config)

    assert exc_info.value.code == 0


def test_minimap2_options_parsing_returns_default_if_none():
    config = {}
    value = analysis_arg_parsing.minimap2_options_parsing(None, config)
    assert isinstance(value, str)
    assert value


def test_minimap2_options_parsing_builds_override_string():
    config = {}
    value = analysis_arg_parsing.minimap2_options_parsing(["k=15", "w=10"], config)
    assert "-k15" in value
    assert "-w10" in value


def test_sample_type_environmental_sets_defaults():
    config = {KEY_SAMPLE_TYPE: "stool"}
    analysis_arg_parsing.sample_type("environmental", config)
    assert config[KEY_SAMPLE_TYPE] == "environmental"

    for key, expected in ENVIRONMENTAL_DEFAULT_DICT.items():
        assert config[key] == expected


def test_analysis_mode_sets_read_length_defaults():
    config = {KEY_ANALYSIS_MODE: "vp1"}
    analysis_arg_parsing.analysis_mode("wg", config)
    assert config[KEY_ANALYSIS_MODE] == "wg"
    assert config[KEY_MIN_READ_LENGTH] is not None
    assert config[KEY_MAX_READ_LENGTH] is not None


def test_haplo_group_parsing_converts_numeric_fields():
    config = {
        KEY_RUN_HAPLOTYPING: False,
        KEY_HAPLOTYPE_SAMPLE_SIZE: "1000",
        KEY_MIN_ALLELE_FREQUENCY: "0.1",
        KEY_MAX_HAPLOTYPES: "5",
        KEY_MIN_HAPLOTYPE_DISTANCE: "2",
        KEY_MIN_HAPLOTYPE_DEPTH: "10",
    }

    analysis_arg_parsing.haplo_group_parsing(
        run_haplotyping=True,
        haplotype_sample_size=2000,
        min_allele_frequency=0.2,
        max_haplotypes=6,
        min_haplotype_distance=3,
        min_haplotype_depth=11,
        config=config,
    )

    assert config[KEY_RUN_HAPLOTYPING] is True
    assert config[KEY_HAPLOTYPE_SAMPLE_SIZE] == 2000
    assert config[KEY_MIN_ALLELE_FREQUENCY] == 0.2
    assert config[KEY_MAX_HAPLOTYPES] == 6
    assert config[KEY_MIN_HAPLOTYPE_DISTANCE] == 3
    assert config[KEY_MIN_HAPLOTYPE_DEPTH] == 11
