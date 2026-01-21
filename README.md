# TieDIE: Tied Diffusion for Subnetwork Discovery

[![Tests](https://github.com/saezlab/tiedie/actions/workflows/ci-testing-unit.yml/badge.svg)](https://github.com/saezlab/tiedie/actions/workflows/ci-testing-unit.yml)
[![Documentation](https://github.com/saezlab/tiedie/actions/workflows/ci-docs.yml/badge.svg)](https://github.com/saezlab/tiedie/actions/workflows/ci-docs.yml)

A network analysis algorithm that finds subnetworks connecting genomic
perturbations to transcriptional changes in large gene interaction networks.

**[Documentation](https://saezlab.github.io/tiedie/)** | **[GitHub](https://github.com/saezlab/tiedie)**

## Authors

Evan O. Paull, Daniel Carlin and Joshua M. Stuart (UC Santa Cruz)

### Additional Contributors

- Srikanth Bezawada (TieDIE Cytoscape Plugin)
- Josh L. Espinoza (Quick kernel loading feature)
- Dana Silverbush (MATLAB kernel generation code updates)
- Denes Turei (modern tooling and packaging)

## Requirements

- Python >= 3.9
- numpy >= 1.20
- scipy >= 1.0
- networkx >= 2.0

## Installation

### Using pip (from GitHub)

```bash
pip install git+https://github.com/saezlab/tiedie.git
```

### Using uv (recommended)

[uv](https://docs.astral.sh/uv/) is a fast Python package manager. To install TieDIE with uv:

```bash
# Install uv if you don't have it
curl -LsSf https://astral.sh/uv/install.sh | sh

# Create a new project with TieDIE
uv init my-project
cd my-project
uv add git+https://github.com/saezlab/tiedie.git

# Or add to an existing project
uv add git+https://github.com/saezlab/tiedie.git
```

### Development installation

```bash
git clone https://github.com/saezlab/tiedie.git
cd tiedie

# Using uv (recommended)
uv venv
source .venv/bin/activate
uv pip install -e ".[tests]"

# Or using pip
python -m venv .venv
source .venv/bin/activate
pip install -e ".[tests]"
```

## Usage

### Command Line Interface

```bash
tiedie -n pathway.sif -u upstream.input -d downstream.input
```

For full options:

```bash
tiedie --help
```

### Python API

```python
from tiedie import ScipyKernel
from tiedie.util import parse_heats, parse_net, normalize_heats

# Parse input files
network = parse_net('pathway.sif')
upstream_heats, upstream_signs = parse_heats('upstream.input')
downstream_heats, downstream_signs = parse_heats('downstream.input')

# Normalize heats
upstream_norm = normalize_heats(upstream_heats)
downstream_norm = normalize_heats(downstream_heats)

# Create diffusion kernel and run diffusion
diffuser = ScipyKernel('pathway.sif')
up_diffused = diffuser.diffuse(upstream_norm, reverse=False)
down_diffused = diffuser.diffuse(downstream_norm, reverse=True)
```

## Input File Formats

### Network file (.sif)

Tab-separated: `source \t interaction \t target`

```
GeneA	-a>	GeneB
GeneB	-t|	GeneC
```

Interaction types:
- `-a>` : activation
- `-t|` : inhibition
- `-component>` : component relationship

### Heat file (.input)

Tab-separated: `gene \t heat \t sign`

```
GeneA	10.5	+
GeneB	8.2	-
```

## Project Structure

```
tiedie/
├── __init__.py      # Public API
├── cli.py           # Command-line interface
├── kernel.py        # Pre-computed kernel diffusion
├── kernel_scipy.py  # On-the-fly kernel generation (scipy)
├── ppr.py           # Personalized PageRank diffusion
├── permute.py       # Permutation testing
├── util.py          # Utility functions
└── ...
tests/               # Test suite
docs/                # Documentation (Tutorial.pdf, FAQ.txt)
```

## Publications

TieDIE was first featured in the 2013 Nature paper "Comprehensive molecular
characterization of clear cell renal cell carcinoma" (TCGA Network). In this
publication, TieDIE analysis connected frequently mutated genes involving the
SWI/SNF chromatin remodelling complex to gene expression changes characteristic
of tumor development and progression. The TieDIE network solution is shown in
Figure 4: http://www.nature.com/nature/journal/v499/n7456/full/nature12222.html

## License

GPL-3.0-or-later

## Contact

Feature requests, comments and requests for clarification should be sent to
<epaull@soe.ucsc.edu>.

## Development

### Running tests

```bash
# Install test dependencies
uv pip install -e ".[tests]"

# Run tests
pytest

# Run tests with coverage
pytest --cov=tiedie
```

### Linting and formatting

The project uses [ruff](https://docs.astral.sh/ruff/) for linting and formatting:

```bash
# Check for issues
ruff check .

# Auto-fix issues
ruff check --fix .

# Format code
ruff format .
```

### Pre-commit hooks

[Pre-commit](https://pre-commit.com/) hooks are configured for automated checks:

```bash
# Install pre-commit hooks
uv pip install -e ".[dev]"
pre-commit install

# Run hooks manually on all files
pre-commit run --all-files
```

### Building documentation

Documentation is built with [MkDocs](https://www.mkdocs.org/):

```bash
# Install docs dependencies
uv pip install -e ".[docs]"

# Serve docs locally
mkdocs serve

# Build static docs
mkdocs build
```

## Version 2.0 Notes

This version represents a major modernization of the TieDIE codebase:

- **Python 3 compatibility**: Ported from Python 2.7 to Python 3.9+
- **Modern packaging**: Converted to a proper Python package with `pyproject.toml`
  and hatchling build backend
- **Test suite**: Added pytest-based test suite with integration tests
- **uv support**: Compatible with modern Python tooling including uv
- **CLI improvements**: Migrated from optparse to argparse

This modernization was done at [Saez Lab](https://saezlab.org/) by
Dénes Türei (<turei.denes@gmail.com>).
