"""Tests for utility functions."""
import doctest

from tiedie import util


def test_doctests() -> None:
    """Run doctests from util module."""
    results = doctest.testmod(util)
    assert results.failed == 0
