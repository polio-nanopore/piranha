import ast
from pathlib import Path
import sys
import types

import pytest


TARGET_PACKAGE_DIRS = [
    "utils",
    "input_parsing",
    "report",
    "analysis",
]


REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))


try:
    import snakemake.common as _snakemake_common

    if not hasattr(_snakemake_common, "DYNAMIC_FILL"):
        _snakemake_common.DYNAMIC_FILL = ""
    if not hasattr(_snakemake_common, "Mode"):
        class _Mode:
            default = "default"
        _snakemake_common.Mode = _Mode
except ModuleNotFoundError:
    snakemake_module = types.ModuleType("snakemake")
    common_module = types.ModuleType("snakemake.common")
    common_module.DYNAMIC_FILL = ""

    class _Mode:
        default = "default"

    common_module.Mode = _Mode
    snakemake_module.common = common_module
    sys.modules["snakemake"] = snakemake_module
    sys.modules["snakemake.common"] = common_module


def _iter_python_files() -> list[Path]:
    repo_root = Path(__file__).resolve().parents[1]
    files: list[Path] = []
    for package_dir in TARGET_PACKAGE_DIRS:
        files.extend(sorted((repo_root / "piranha" / package_dir).glob("*.py")))
    return [path for path in files if path.name != "__init__.py"]


def discover_target_functions() -> list[tuple[str, str, str]]:
    discovered: list[tuple[str, str, str]] = []

    for path in _iter_python_files():
        module_name = ".".join(path.with_suffix("").parts[-3:])
        tree = ast.parse(path.read_text(), filename=str(path))

        for node in tree.body:
            if isinstance(node, ast.FunctionDef):
                discovered.append((module_name, node.name, str(path)))

    return discovered


@pytest.fixture(scope="session")
def all_target_functions() -> list[tuple[str, str, str]]:
    return discover_target_functions()
