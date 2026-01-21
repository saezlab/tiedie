"""Tests for utility functions."""

import pytest

from tiedie.util import (
    parse_net,
    parse_heats,
    normalize_heats,
    get_network_nodes,
)


class TestUtilFunctions:
    """Tests for utility parsing and normalization functions."""

    def test_parse_heats(self, tmp_path: pytest.fixture) -> None:
        """Test parse_heats function with simple input."""
        heat_file = tmp_path / 'test.heats'
        heat_file.write_text('A\t10\t+\nB\t20\t-\n')

        heats, signs = parse_heats(str(heat_file))

        assert heats == {'A': 10.0, 'B': 20.0}
        assert signs == {'A': '+', 'B': '-'}

    def test_parse_net(self, tmp_path: pytest.fixture) -> None:
        """Test parse_net function with simple SIF input."""
        net_file = tmp_path / 'test.sif'
        net_file.write_text('A\t-a>\tB\nB\t-t|\tC\n')

        network = parse_net(str(net_file))

        assert 'A' in network
        assert 'B' in network
        # A should have edge to B
        assert any(target == 'B' for _, target in network['A'])

    def test_get_network_nodes(self, tmp_path: pytest.fixture) -> None:
        """Test get_network_nodes returns all nodes."""
        net_file = tmp_path / 'test.sif'
        net_file.write_text('A\t-a>\tB\nB\t-t|\tC\n')

        network = parse_net(str(net_file))
        nodes = get_network_nodes(network)

        assert 'A' in nodes
        assert 'B' in nodes
        assert 'C' in nodes

    def test_normalize_heats(self) -> None:
        """Test normalize_heats normalization."""
        heats = {'A': 100.0, 'B': 100.0}

        normalized = normalize_heats(heats)

        # Sum of absolute values should equal 1000 (FACTOR in normalize_heats)
        total = sum(abs(v) for v in normalized.values())
        assert total == pytest.approx(1000.0)
