"""Pytest configuration and fixtures."""

from pathlib import Path

import pytest

TEST_DIR = Path(__file__).parent / 'test_files'


@pytest.fixture
def test_dir() -> Path:
    """Return the path to the test files directory."""
    return TEST_DIR


@pytest.fixture
def gbm_pathway(test_dir: Path) -> Path:
    """Return path to GBM pathway file."""
    return test_dir / 'gbm_pathway.sif'


@pytest.fixture
def gbm_upstream(test_dir: Path) -> Path:
    """Return path to GBM upstream heats file."""
    return test_dir / 'gbm_upstream.input'


@pytest.fixture
def gbm_downstream(test_dir: Path) -> Path:
    """Return path to GBM downstream heats file."""
    return test_dir / 'gbm_downstream.input'
