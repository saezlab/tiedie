"""Integration tests for TieDIE using GBM example data."""
from pathlib import Path

from tiedie import SciPYKernel
from tiedie.util import getNetworkNodes, normalizeHeats, parseHeats, parseNet

TEST_DIR = Path(__file__).parent / 'test_files'
GBM_PATHWAY = TEST_DIR / 'gbm_pathway.sif'
GBM_UPSTREAM = TEST_DIR / 'gbm_upstream.input'
GBM_DOWNSTREAM = TEST_DIR / 'gbm_downstream.input'


class TestGBMExample:
    """Integration tests using the GBM (Glioblastoma) example data."""

    def test_parse_network(self) -> None:
        """Test parsing the GBM pathway network."""
        network = parseNet(str(GBM_PATHWAY))

        assert len(network) > 0
        # Check that known nodes are present
        nodes = getNetworkNodes(network)
        assert 'PIK3CA' in nodes
        assert 'PTEN' in nodes
        assert 'AKT1' in nodes

    def test_parse_heats(self) -> None:
        """Test parsing upstream and downstream heat inputs."""
        upstream_heats, upstream_signs = parseHeats(str(GBM_UPSTREAM))
        downstream_heats, downstream_signs = parseHeats(str(GBM_DOWNSTREAM))

        assert len(upstream_heats) > 0
        assert len(downstream_heats) > 0
        assert 'PIK3R1' in upstream_heats
        assert 'PTEN' in upstream_heats
        assert 'G2/M_arrest' in downstream_heats

    def test_diffuse_heats(self) -> None:
        """Test heat diffusion on the GBM network."""
        upstream_heats, _ = parseHeats(str(GBM_UPSTREAM))

        # Normalize heats
        normalized = normalizeHeats(upstream_heats)

        # Create kernel and diffuse
        diffuser = SciPYKernel(str(GBM_PATHWAY))
        diffused = diffuser.diffuse(normalized, reverse=False)

        # Check that diffusion produces non-zero heats
        assert len(diffused) > 0
        assert any(v > 0 for v in diffused.values())

        # Nodes with initial heat should have higher diffused values
        for node in normalized:
            if normalized[node] > 0 and node in diffused:
                assert diffused[node] > 0

    def test_bidirectional_diffusion(self) -> None:
        """Test that bidirectional diffusion produces connecting nodes."""
        upstream_heats, _ = parseHeats(str(GBM_UPSTREAM))
        downstream_heats, _ = parseHeats(str(GBM_DOWNSTREAM))

        # Normalize heats
        upstream_norm = normalizeHeats(upstream_heats)
        downstream_norm = normalizeHeats(downstream_heats)

        # Create kernel
        diffuser = SciPYKernel(str(GBM_PATHWAY))

        # Diffuse from upstream (forward)
        up_diffused = diffuser.diffuse(upstream_norm, reverse=False)

        # Diffuse from downstream (reverse)
        down_diffused = diffuser.diffuse(downstream_norm, reverse=True)

        # Check that both directions produce results
        assert len(up_diffused) > 0
        assert len(down_diffused) > 0

        # Find nodes with heat from both directions (potential linkers)
        linker_candidates = set(up_diffused.keys()) & set(down_diffused.keys())
        assert len(linker_candidates) > 0
