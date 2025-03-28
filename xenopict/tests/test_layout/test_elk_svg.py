"""Test ELK SVG layout functionality."""

from lxml import etree

from xenopict.layout.elk import layout_to_svg


def _parse_svg(svg_str: str) -> etree._Element:
    """Parse SVG string into an XML element tree."""
    parser = etree.XMLParser(remove_blank_text=True)
    return etree.fromstring(svg_str.encode(), parser=parser)


def test_simple_svg_layout():
    """Test basic graph layout to SVG conversion."""
    graph = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30},
            {"id": "n2", "width": 30, "height": 30},
        ],
        "edges": [{"id": "e1", "sources": ["n1"], "targets": ["n2"]}],
    }

    svg = layout_to_svg(graph)
    root = _parse_svg(svg)

    # Basic SVG validation
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert root.get("version") == "1.1"

    # Find elements by class
    nodes = root.findall(".//*[@class='node']")
    edges = root.findall(".//*[@class='edge']")
    assert len(nodes) == 2  # Should have two nodes
    assert len(edges) == 1  # Should have one edge

    # Check element IDs
    node_ids = {node.get("id") for node in nodes}
    edge_ids = {edge.get("id") for edge in edges}
    assert node_ids == {"n1", "n2"}
    assert edge_ids == {"e1"}


def test_svg_layout_with_options():
    """Test SVG layout with custom options."""
    graph = {
        "id": "root",
        "children": [
            {"id": "n1", "width": 30, "height": 30},
            {"id": "n2", "width": 30, "height": 30},
        ],
        "edges": [{"id": "e1", "sources": ["n1"], "targets": ["n2"]}],
    }

    options = {"elk.algorithm": "stress", "elk.spacing.nodeNode": "50"}
    svg = layout_to_svg(graph, options)
    root = _parse_svg(svg)

    # Basic SVG validation
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert root.get("version") == "1.1"

    # Find elements by class
    nodes = root.findall(".//*[@class='node']")
    edges = root.findall(".//*[@class='edge']")
    assert len(nodes) == 2  # Should have two nodes
    assert len(edges) == 1  # Should have one edge

    # Check element IDs
    node_ids = {node.get("id") for node in nodes}
    edge_ids = {edge.get("id") for edge in edges}
    assert node_ids == {"n1", "n2"}
    assert edge_ids == {"e1"}


def test_complex_svg_layout():
    """Test SVG layout with a more complex graph structure."""
    graph = {
        "id": "root",
        "children": [
            {
                "id": "compound1",
                "children": [
                    {"id": "n1", "width": 30, "height": 30},
                    {"id": "n2", "width": 30, "height": 30},
                ],
                "edges": [{"id": "e1", "sources": ["n1"], "targets": ["n2"]}],
            },
            {
                "id": "compound2",
                "children": [
                    {"id": "n3", "width": 30, "height": 30},
                    {"id": "n4", "width": 30, "height": 30},
                ],
                "edges": [{"id": "e2", "sources": ["n3"], "targets": ["n4"]}],
            },
        ],
        "edges": [
            {"id": "e3", "sources": ["n2"], "targets": ["n3"]},
        ],
    }

    svg = layout_to_svg(graph)
    root = _parse_svg(svg)

    # Basic SVG validation
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert root.get("version") == "1.1"

    # Find all elements by ID
    elements = root.findall(".//*[@id]")
    element_ids = {elem.get("id") for elem in elements}

    # Check for compound nodes
    assert "compound1" in element_ids
    assert "compound2" in element_ids

    # Check for all nodes
    for node_id in ["n1", "n2", "n3", "n4"]:
        assert node_id in element_ids

    # Check for all edges
    for edge_id in ["e1", "e2", "e3"]:
        assert edge_id in element_ids


def test_svg_with_empty_graph():
    """Test SVG generation with an empty graph."""
    graph = {
        "id": "root",
        "children": [],
        "edges": [],
    }

    svg = layout_to_svg(graph)
    root = _parse_svg(svg)

    # Basic SVG validation
    assert root.tag == "{http://www.w3.org/2000/svg}svg"
    assert root.get("version") == "1.1"

    # Find root group element
    root_group = root.findall(".//*[@id='root']")
    assert len(root_group) == 1  # Should have exactly one root group element
