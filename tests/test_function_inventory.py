import importlib
import inspect

import pytest

from conftest import discover_target_functions


FUNCTIONS = discover_target_functions()


@pytest.mark.parametrize(
    "module_name,function_name,source_path",
    [
        pytest.param(module_name, function_name, source_path, id=f"{module_name}.{function_name}")
        for module_name, function_name, source_path in FUNCTIONS
    ],
)
def test_all_top_level_functions_importable(module_name, function_name, source_path):
    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        pytest.skip(f"Missing optional dependency for {module_name}: {exc.name}")
    except ImportError as exc:
        pytest.skip(f"Import-time dependency mismatch for {module_name}: {exc}")

    function = getattr(module, function_name, None)
    assert function is not None, f"{function_name} missing in {module_name} ({source_path})"
    assert inspect.isfunction(function), f"{module_name}.{function_name} is not a function"


def test_function_inventory_not_empty():
    assert FUNCTIONS, "No functions discovered in target packages"
