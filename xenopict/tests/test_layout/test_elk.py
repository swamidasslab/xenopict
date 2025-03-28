"""Test ELK layout functionality."""

from xenopict.layout.elk import get_layout_algorithms, get_layout_options, layout


def test_simple_layout():
    """Test basic graph layout with ELK."""
    graph = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30},
            {"id": "n2", "width": 30, "height": 30},
        ],
        "edges": [{"id": "e1", "sources": ["n1"], "targets": ["n2"]}],
    }

    result = layout(graph)

    assert isinstance(result, dict)
    assert "children" in result
    assert len(result["children"]) == 2
    assert isinstance(result["children"][0]["x"], (int, float))
    assert isinstance(result["children"][0]["y"], (int, float))


def test_layout_with_options():
    """Test layout with custom options."""
    graph = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30},
            {"id": "n2", "width": 30, "height": 30},
        ],
        "edges": [{"id": "e1", "sources": ["n1"], "targets": ["n2"]}],
    }

    options = {"elk.algorithm": "stress", "elk.spacing.nodeNode": "50"}

    result = layout(graph, options)

    assert isinstance(result, dict)
    assert "children" in result
    assert len(result["children"]) == 2
    assert isinstance(result["children"][0]["x"], (int, float))
    assert isinstance(result["children"][0]["y"], (int, float))


def test_get_layout_options():
    """Test retrieving available layout options."""
    options = get_layout_options()
    assert isinstance(options, list)
    assert len(options) > 0
    assert all(isinstance(opt, dict) for opt in options)
    assert all("id" in opt for opt in options)


def test_get_layout_algorithms():
    """Test retrieving available layout algorithms."""
    algorithms = get_layout_algorithms()
    assert isinstance(algorithms, list)
    assert len(algorithms) > 0
    assert all(isinstance(algo, dict) for algo in algorithms)
    assert all("id" in algo for algo in algorithms)
    # Check that common algorithms are available
    common_algorithms = {
        "org.eclipse.elk.layered",
        "org.eclipse.elk.force",
        "org.eclipse.elk.radial",
    }
    found_algorithms = set(algo["id"] for algo in algorithms)
    assert common_algorithms.intersection(found_algorithms)
