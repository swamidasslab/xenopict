"""Tests for the ELK layout integration."""

import pytest
from xenopict.layout.elk import layout, get_layout_options, get_layout_algorithms

def test_simple_layout():
    """Test basic graph layout with ELK."""
    graph = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30},
            {"id": "n2", "width": 30, "height": 30}
        ],
        "edges": [
            {"id": "e1", "sources": ["n1"], "targets": ["n2"]}
        ]
    }
    
    result = layout(graph)
    
    # Check that layout was applied
    assert "children" in result
    assert len(result["children"]) == 2
    
    # Check that positions were assigned
    for node in result["children"]:
        assert "x" in node
        assert "y" in node
        assert isinstance(node["x"], (int, float))
        assert isinstance(node["y"], (int, float))

def test_layout_with_options():
    """Test layout with custom options."""
    graph = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30},
            {"id": "n2", "width": 30, "height": 30}
        ],
        "edges": [
            {"id": "e1", "sources": ["n1"], "targets": ["n2"]}
        ]
    }
    
    options = {
        "elk.algorithm": "stress",
        "elk.spacing.nodeNode": "50"
    }
    
    result = layout(graph, options)
    
    # Check that layout was applied
    assert "children" in result
    assert len(result["children"]) == 2
    
    # Check that positions were assigned
    for node in result["children"]:
        assert "x" in node
        assert "y" in node
        assert isinstance(node["x"], (int, float))
        assert isinstance(node["y"], (int, float))

def test_get_layout_options():
    """Test retrieving available layout options."""
    options = get_layout_options()
    assert isinstance(options, dict)
    assert len(options) > 0

def test_get_layout_algorithms():
    """Test retrieving available layout algorithms."""
    algorithms = get_layout_algorithms()
    assert isinstance(algorithms, dict)
    assert len(algorithms) > 0
    # Check that common algorithms are available
    common_algorithms = {"layered", "stress", "force"}
    found_algorithms = set(algorithms.keys())
    assert common_algorithms.intersection(found_algorithms) 