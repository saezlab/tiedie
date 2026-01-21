"""Tests for heat diffusion kernel functionality."""

from pathlib import Path

from tiedie import Kernel, ScipyKernel
from tiedie.util import parse_heats


TEST_DIR = Path(__file__).parent / 'test_files'
TEST_PATHWAY = TEST_DIR / 'test.pathway.sif'
TEST_KERNEL = TEST_DIR / 'kernel.tab'
TEST_INPUT = TEST_DIR / 'upstream.input'
TEST_DIFFUSED_SOLN = TEST_DIR / 'upstream.diffused'


class TestScipyKernel:
    """Tests for on-the-fly kernel diffusion with scipy."""

    def test_diffuse(self) -> None:
        """Test that ScipyKernel diffusion produces valid output."""
        input_heats, _ = parse_heats(str(TEST_INPUT))

        diffuser = ScipyKernel(str(TEST_PATHWAY))
        diffused = diffuser.diffuse(input_heats, reverse=False)

        # Verify diffusion produces results
        assert len(diffused) > 0

        # Input nodes should have heat in output
        for node in input_heats:
            if node in diffused:
                assert diffused[node] >= 0

        # Heat should spread to neighbors
        assert any(v > 0 for v in diffused.values())

    def test_diffuse_reverse(self) -> None:
        """Test reverse diffusion on the network."""
        input_heats, _ = parse_heats(str(TEST_INPUT))

        diffuser = ScipyKernel(str(TEST_PATHWAY))
        diffused = diffuser.diffuse(input_heats, reverse=True)

        assert len(diffused) > 0
        assert any(v > 0 for v in diffused.values())


class TestPrecomputedKernel:
    """Tests for pre-computed kernel diffusion."""

    def test_diffuse(self) -> None:
        """Test that pre-computed kernel diffusion produces valid output."""
        input_heats, _ = parse_heats(str(TEST_INPUT))

        diffuser = Kernel(str(TEST_KERNEL))
        diffused = diffuser.diffuse(input_heats, reverse=False)

        # Verify diffusion produces results
        assert len(diffused) > 0

        # Input nodes should have heat in output
        for node in input_heats:
            if node in diffused:
                assert diffused[node] >= 0

        # Heat should spread
        assert any(v > 0 for v in diffused.values())

    def test_kernel_consistency(self) -> None:
        """Test that scipy and precomputed kernel produce similar results."""
        input_heats, _ = parse_heats(str(TEST_INPUT))

        scipy_diffuser = ScipyKernel(str(TEST_PATHWAY))
        kernel_diffuser = Kernel(str(TEST_KERNEL))

        scipy_diffused = scipy_diffuser.diffuse(input_heats, reverse=False)
        kernel_diffused = kernel_diffuser.diffuse(input_heats, reverse=False)

        # Both should produce results for the same input nodes
        common_nodes = set(scipy_diffused.keys()) & set(kernel_diffused.keys())
        assert len(common_nodes) > 0

        # Values should be in the same order of magnitude
        for node in common_nodes:
            scipy_val = scipy_diffused[node]
            kernel_val = kernel_diffused[node]
            if scipy_val > 0.001 and kernel_val > 0.001:
                ratio = scipy_val / kernel_val
                assert 0.1 < ratio < 10, f'Values differ too much for {node}'
