import pytest

from piranha.utils import dependency_checks, log_colours, misc


def test_log_colour_wrappers_add_escape_codes():
    assert log_colours.red("x").startswith(log_colours.RED)
    assert log_colours.green("x").startswith(log_colours.GREEN)
    assert log_colours.cyan("x").startswith(log_colours.CYAN)
    assert log_colours.yellow("x").startswith(log_colours.YELLOW)


def test_add_arg_to_config_updates_for_truthy_and_zero():
    config = {"threads": 1}
    misc.add_arg_to_config("threads", 4, config)
    misc.add_arg_to_config("min_reads", 0, config)

    assert config["threads"] == 4
    assert config["min_reads"] == 0


def test_check_date_format_accepts_iso_date():
    misc.check_date_format("2026-03-26", 1, "date")


def test_check_date_format_rejects_bad_date():
    with pytest.raises(SystemExit):
        misc.check_date_format("26-03-2026", 1, "date")


def test_dependency_check_module_appends_missing_module():
    missing = []
    dependency_checks.check_module("definitely_not_a_real_module_123", missing)
    assert "definitely_not_a_real_module_123" in missing


def test_check_dependencies_exits_on_missing_dependency(monkeypatch):
    monkeypatch.setattr(dependency_checks, "which", lambda dep: False)

    with pytest.raises(SystemExit):
        dependency_checks.check_dependencies(["nonexistent_binary"], [])
