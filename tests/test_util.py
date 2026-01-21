"""Tests for utility functions."""

import pytest

from tiedie.util import parseNet, parseHeats, normalizeHeats, getNetworkNodes


class TestUtilFunctions:
    """Tests for utility parsing and normalization functions."""

    def test_parse_heats(self, tmp_path: pytest.fixture) -> None:
        """Test parseHeats function with simple input."""
        heat_file = tmp_path / 'test.heats'
        heat_file.write_text('A\t10\t+\nB\t20\t-\n')

        heats, signs = parseHeats(str(heat_file))

        assert heats == {'A': 10.0, 'B': 20.0}
        assert signs == {'A': '+', 'B': '-'}

    def test_parse_net(self, tmp_path: pytest.fixture) -> None:
        """Test parseNet function with simple SIF input."""
        net_file = tmp_path / 'test.sif'
        net_file.write_text('A\t-a>\tB\nB\t-t|\tC\n')

        network = parseNet(str(net_file))

        assert 'A' in network
        assert 'B' in network
        # A should have edge to B
        assert any(target == 'B' for _, target in network['A'])

    def test_get_network_nodes(self, tmp_path: pytest.fixture) -> None:
        """Test getNetworkNodes returns all nodes."""
        net_file = tmp_path / 'test.sif'
        net_file.write_text('A\t-a>\tB\nB\t-t|\tC\n')

        network = parseNet(str(net_file))
        nodes = getNetworkNodes(network)

        assert 'A' in nodes
        assert 'B' in nodes
        assert 'C' in nodes

    def test_normalize_heats(self) -> None:
        """Test normalizeHeats normalization."""
        heats = {'A': 100.0, 'B': 100.0}

        normalized = normalizeHeats(heats)

        # Sum of absolute values should equal 1000 (FACTOR in normalizeHeats)
        total = sum(abs(v) for v in normalized.values())
        assert total == pytest.approx(1000.0)
