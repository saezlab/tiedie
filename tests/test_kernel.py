"""Tests for heat diffusion kernel functionality."""
from pathlib import Path

import pytest

from tiedie import Kernel, SciPYKernel
from tiedie.util import parseHeats

TEST_DIR = Path(__file__).parent / 'test_files'
TEST_PATHWAY = TEST_DIR / 'test.pathway.sif'
TEST_KERNEL = TEST_DIR / 'kernel.tab'
TEST_INPUT = TEST_DIR / 'upstream.input'
TEST_DIFFUSED_SOLN = TEST_DIR / 'upstream.diffused'


class TestSciPYKernel:
    """Tests for on-the-fly kernel diffusion with scipy."""

    def test_diffuse(self) -> None:
        """Test that SciPYKernel diffusion matches expected output."""
        correct_heats, _ = parseHeats(str(TEST_DIFFUSED_SOLN))
        input_heats, _ = parseHeats(str(TEST_INPUT))

        diffuser = SciPYKernel(str(TEST_PATHWAY))
        diffused = diffuser.diffuse(input_heats, reverse=False)

        for key, val in diffused.items():
            assert val == pytest.approx(correct_heats[key], abs=1e-10)


class TestPrecomputedKernel:
    """Tests for pre-computed kernel diffusion."""

    def test_diffuse(self) -> None:
        """Test that pre-computed kernel diffusion matches expected output."""
        correct_heats, _ = parseHeats(str(TEST_DIFFUSED_SOLN))
        input_heats, _ = parseHeats(str(TEST_INPUT))

        diffuser = Kernel(str(TEST_KERNEL))
        diffused = diffuser.diffuse(input_heats, reverse=False)

        for key, val in diffused.items():
            # Lower precision for higher values due to numerical differences
            if val < 1:
                assert val == pytest.approx(correct_heats[key], abs=1e-5)
            else:
                assert val == pytest.approx(correct_heats[key], abs=1e-3)
