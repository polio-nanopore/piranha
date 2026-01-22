# Development Setup Guide

## Project Structure

Piranha uses modern Python packaging standards with configuration centralized in `pyproject.toml`.

### Key Files
- **pyproject.toml** - Main project configuration (build, dependencies, tool settings)
- **setup.py** - Minimal backwards-compatibility shim
- **setup.cfg** - Fallback configuration for older tools
- **tox.ini** - Testing and linting configuration
- **MANIFEST.in** - Specification of non-Python files to include in distributions

## Installation for Development

### 1. Clone and Setup Environment

```bash
git clone https://github.com/polio-nanopore/piranha.git
cd piranha
```

### 2. Create Development Environment (with conda/mamba)

```bash
mamba env create -f environment.yml
conda activate piranha
pip install -e ".[dev]"
```

Or with pip only:

```bash
python -m venv venv
source venv/bin/activate
pip install -e ".[dev]"
```

### 3. Setup Pre-commit Hooks (Recommended)

```bash
pip install pre-commit
pre-commit install
```

This will automatically run code quality checks before each commit.

## Code Quality Tools

### Black (Code Formatting)
```bash
black piranha
```

### isort (Import Sorting)
```bash
isort piranha
```

### Ruff (Fast Linting)
```bash
ruff check piranha
ruff check --fix piranha  # Auto-fix issues
```

### Flake8 (Style Checking)
```bash
flake8 piranha
```

### MyPy (Type Checking)
```bash
mypy piranha
```

### Running All Checks
```bash
tox -e lint    # Run all linting checks
tox -e type    # Run type checking
tox -e format  # Auto-format code
```

## Testing

### Run Tests
```bash
pytest
pytest -v  # Verbose output
pytest --cov=piranha  # With coverage
```

### Test Specific Module
```bash
pytest piranha/test/test_consensus_functions.py
```

### Generate Coverage Report
```bash
pytest --cov=piranha --cov-report=html
open htmlcov/index.html
```

### Test Across Python Versions (requires tox)
```bash
tox
```

## Type Hints

We recommend gradually adding type hints to the codebase. Here's how:

```python
from typing import List, Dict, Optional, Tuple
from pathlib import Path

def process_sequences(sequences: List[str], output_dir: Path) -> Dict[str, str]:
    """Process viral sequences.
    
    Args:
        sequences: List of FASTA sequences
        output_dir: Directory to write outputs
        
    Returns:
        Dictionary mapping sequence IDs to processed sequences
    """
    results: Dict[str, str] = {}
    # implementation
    return results
```

## Documentation

### Building Docs (when Sphinx is added)
```bash
pip install ".[docs]"
cd docs
make html
open _build/html/index.html
```

### Docstring Style

Use Google-style docstrings:

```python
def find_variants(reference_seq: str, query_seq: str) -> List[str]:
    """Identify variants between reference and query sequences.
    
    Args:
        reference_seq: Reference sequence
        query_seq: Query sequence to compare
        
    Returns:
        List of variant annotations in format "position:refquery"
        
    Raises:
        ValueError: If sequences have different lengths
    """
```

## Publishing

### Build Distribution
```bash
pip install build
python -m build
```

### Upload to PyPI (requires credentials)
```bash
pip install twine
twine upload dist/*
```

## Dependency Management

### Adding Dependencies

1. Add to `pyproject.toml` under `dependencies` (required) or `[project.optional-dependencies]` (optional)
2. Ensure conda environment users update with: `mamba env update -f environment.yml`
3. Update version constraints if needed

### Version Constraints in pyproject.toml

- `package==1.2.3` - Exact version
- `package>=1.2.3` - Minimum version
- `package~=1.2` - Compatible release (~= 1.2 means >= 1.2, < 2.0)
- `package>=1.2,<2.0` - Version range

## Troubleshooting

### Import Errors After Install
```bash
pip install -e . --force-reinstall
```

### Type Checking Issues
```bash
mypy piranha --show-error-codes
```

### Test Discovery Issues
Ensure test files follow pattern: `test_*.py` and are in `piranha/test/`

## Useful Resources

- [Python Packaging Guide](https://packaging.python.org/)
- [PEP 517 - Build Backend](https://peps.python.org/pep-0517/)
- [PEP 518 - pyproject.toml](https://peps.python.org/pep-0518/)
- [Black Documentation](https://black.readthedocs.io/)
- [MyPy Documentation](https://mypy.readthedocs.io/)
