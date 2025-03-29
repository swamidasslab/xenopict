"""Test ELK layout functionality."""

from xenopict.layout.elk import (
    EdgeRouting,
    LayoutAlgorithm,
    LayoutConfig,
    LayoutDirection,
    LayoutPadding,
    LayoutSpacing,
    NodeSizeConfig,
    _build_layout_options,
    elk_to_layout,
    elk_to_svg,
    get_layout_algorithms,
    get_layout_options,
    layout_to_svg,
)


def test_elk_to_layout():
    """Test converting ELK graph to layout."""
    graph = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30, "x": 0, "y": 0},
            {"id": "n2", "width": 30, "height": 30, "x": 100, "y": 0},
        ],
        "edges": [
            {
                "id": "e1",
                "sources": ["n1"],
                "targets": ["n2"],
                "sections": [
                    {
                        "startPoint": {"x": 30, "y": 15},
                        "endPoint": {"x": 100, "y": 15},
                    }
                ],
            }
        ],
    }
    layout = elk_to_layout(graph)
    assert layout is not None


def test_node_size_defaults():
    """Test that node sizes are handled correctly with defaults and overrides."""
    # Test graph with a mix of default and explicit sizes
    graph = {
        "id": "root",
        "children": [
            {"id": "n1"},  # Should get default size
            {"id": "n2", "width": 150, "height": 75},  # Explicit size
            {"id": "n3", "width": 30, "height": 30},  # Below minimum size
        ],
        "edges": [
            {"id": "e1", "sources": ["n1"], "targets": ["n2"]},
            {"id": "e2", "sources": ["n2"], "targets": ["n3"]},
        ],
    }

    # Test with default options
    result = elk_to_layout(graph)
    print("Default options result:", result)

    # Check node with default size (n1)
    n1 = next(n for n in result["children"] if n["id"] == "n1")
    assert n1["width"] == 100, "Default width should be 100"
    assert n1["height"] == 100, "Default height should be 100"

    # Check node with explicit size (n2)
    n2 = next(n for n in result["children"] if n["id"] == "n2")
    assert n2["width"] == 150, "Explicit width should be preserved"
    assert n2["height"] == 75, "Explicit height should be preserved"

    # Check node with size below minimum (n3)
    n3 = next(n for n in result["children"] if n["id"] == "n3")
    assert n3["width"] == 50, "Width should be adjusted to minimum"
    assert n3["height"] == 50, "Height should be adjusted to minimum"

    # Test with custom node sizes in the graph
    graph_with_custom = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 200, "height": 150},  # Custom size
            {"id": "n2", "width": 150, "height": 75},  # Explicit size
            {"id": "n3", "width": 80, "height": 60},  # Above minimum size
        ],
        "edges": [
            {"id": "e1", "sources": ["n1"], "targets": ["n2"]},
            {"id": "e2", "sources": ["n2"], "targets": ["n3"]},
        ],
    }
    result = elk_to_layout(graph_with_custom)
    print("Custom sizes result:", result)

    # Check node with custom size (n1)
    n1 = next(n for n in result["children"] if n["id"] == "n1")
    assert n1["width"] == 200, "Custom width should be preserved"
    assert n1["height"] == 150, "Custom height should be preserved"

    # Check node with explicit size (n2)
    n2 = next(n for n in result["children"] if n["id"] == "n2")
    assert n2["width"] == 150, "Explicit width should still be preserved"
    assert n2["height"] == 75, "Explicit height should still be preserved"

    # Check node with size above minimum (n3)
    n3 = next(n for n in result["children"] if n["id"] == "n3")
    assert n3["width"] == 80, "Width above minimum should be preserved"
    assert n3["height"] == 60, "Height above minimum should be preserved"


def test_elk_to_svg():
    """Test converting ELK graph to SVG."""
    graph = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30, "x": 0, "y": 0},
            {"id": "n2", "width": 30, "height": 30, "x": 100, "y": 0},
        ],
        "edges": [
            {
                "id": "e1",
                "sources": ["n1"],
                "targets": ["n2"],
                "sections": [
                    {
                        "startPoint": {"x": 30, "y": 15},
                        "endPoint": {"x": 100, "y": 15},
                    }
                ],
            }
        ],
    }
    svg = elk_to_svg(graph)
    assert svg is not None
    assert isinstance(svg, str)
    assert svg.startswith("<?xml")
    assert svg.endswith("</svg>")


def test_layout_to_svg():
    """Test converting layout to SVG."""
    layout = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30, "x": 0, "y": 0},
            {"id": "n2", "width": 30, "height": 30, "x": 100, "y": 0},
        ],
        "edges": [
            {
                "id": "e1",
                "sources": ["n1"],
                "targets": ["n2"],
                "sections": [
                    {
                        "startPoint": {"x": 30, "y": 15},
                        "endPoint": {"x": 100, "y": 15},
                    }
                ],
            }
        ],
    }
    svg = layout_to_svg(layout)
    assert svg is not None
    assert isinstance(svg, str)
    assert svg.startswith("<?xml")
    assert svg.endswith("</svg>")


# def test_layout_with_options():
#     """Test layout with custom options."""
#     graph = {
#         "id": "root",
#         "children": [
#             {"id": "n1", "width": 30, "height": 30},
#             {"id": "n2", "width": 30, "height": 30},
#         ],
#         "edges": [{"id": "e1", "sources": ["n1"], "targets": ["n2"]}],
#     }

#     options = {"elk.algorithm": "stress", "elk.spacing.nodeNode": "50"}

#     result = layout(graph, options)

#     assert isinstance(result, dict)
#     assert "children" in result
#     assert len(result["children"]) == 2
#     assert isinstance(result["children"][0]["x"], (int, float))
#     assert isinstance(result["children"][0]["y"], (int, float))


def test_get_layout_options():
    """Test getting layout options."""
    assert get_layout_options() is not None


def test_get_layout_algorithms():
    """Test getting layout algorithms."""
    assert get_layout_algorithms() is not None


def test_layout_config_with_options():
    """Test that LayoutConfig works correctly with actual ELK options."""
    # Create a config with custom values
    config = LayoutConfig(
        algorithm=LayoutAlgorithm.STRESS,
        direction=LayoutDirection.DOWN,
        edge_routing=EdgeRouting.POLYLINE,
        aspect_ratio=2.0,
        node_size=NodeSizeConfig(
            default_width=200, default_height=150, min_width=75, min_height=60
        ),
        spacing=LayoutSpacing(node_node=30, edge_edge=15, edge_node=20),
        padding=LayoutPadding(top=30, left=30, bottom=30, right=30),
    )

    # Convert config to ELK options
    elk_options = _build_layout_options(config)

    # Test that generated options match ELK's format
    assert elk_options["elk.algorithm"] == "stress"
    assert elk_options["elk.direction"] == "DOWN"
    assert elk_options["elk.edgeRouting"] == "POLYLINE"
    assert float(elk_options["elk.aspectRatio"]) == 2.0
    assert elk_options["elk.spacing.nodeNode"] == "30"
    assert elk_options["elk.spacing.edgeEdge"] == "15"
    assert elk_options["elk.spacing.edgeNode"] == "20"
    assert elk_options["elk.padding"] == "top=30,left=30,bottom=30,right=30"

    # Test with a graph to ensure options are applied
    graph = {
        "id": "root",
        "children": [
            {"id": "n1"},  # Should get default size
            {"id": "n2", "width": 50, "height": 40},  # Should be adjusted to minimum
        ],
        "edges": [{"id": "e1", "sources": ["n1"], "targets": ["n2"]}],
    }

    result = elk_to_layout(graph, config)

    # Check that node sizes were applied correctly
    n1 = next(n for n in result["children"] if n["id"] == "n1")
    assert n1["width"] == 200, "Default width should be applied"
    assert n1["height"] == 150, "Default height should be applied"

    n2 = next(n for n in result["children"] if n["id"] == "n2")
    assert n2["width"] == 75, "Width should be adjusted to minimum"
    assert n2["height"] == 60, "Height should be adjusted to minimum"
