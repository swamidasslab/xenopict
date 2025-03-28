"""Test SVG integration with ELK layout."""

from lxml import etree

from xenopict.layout.svg_elk import (
    SVG_NS,
    create_elk_graph,
    get_svg_size,
    layout_with_svgs,
    render_layout_svg,
)

# Create a parser that handles SVG
parser = etree.XMLParser(remove_blank_text=True)


def test_get_svg_size():
    """Test SVG size extraction."""
    # Test with explicit width/height
    svg = etree.fromstring('<svg width="100" height="50"></svg>', parser=parser)
    width, height = get_svg_size(svg)
    assert width == 100
    assert height == 50

    # Test with viewBox
    svg = etree.fromstring('<svg viewBox="0 0 200 100"></svg>', parser=parser)
    width, height = get_svg_size(svg)
    assert width == 200
    assert height == 100


def test_create_elk_graph():
    """Test ELK graph creation from SVGs."""
    # Create test SVGs
    nodes = {
        "n1": etree.fromstring('<svg width="100" height="50"></svg>', parser=parser),
        "n2": etree.fromstring('<svg width="120" height="60"></svg>', parser=parser),
    }
    edges = {"e1": ("n1", "n2")}

    graph = create_elk_graph(nodes, edges)

    # Verify graph structure
    assert graph["id"] == "root"
    assert len(graph["children"]) == 2
    assert len(graph["edges"]) == 1

    # Verify node properties
    node1 = next(n for n in graph["children"] if n["id"] == "n1")
    assert node1["width"] == 100
    assert node1["height"] == 50

    # Verify edge properties
    edge = graph["edges"][0]
    assert edge["id"] == "e1"
    assert edge["sources"] == ["n1"]
    assert edge["targets"] == ["n2"]


def test_render_layout_svg():
    """Test rendering layout with SVG nodes."""
    # Create test SVGs with different shapes
    nodes = {
        "n1": etree.fromstring(
            '<svg width="100" height="50">'
            '<rect width="80" height="30" fill="blue"/>'
            "</svg>",
            parser=parser,
        ),
        "n2": etree.fromstring(
            '<svg width="120" height="60">'
            '<circle r="25" cx="60" cy="30" fill="red"/>'
            "</svg>",
            parser=parser,
        ),
    }

    # Test with edges and bend points
    layout_result = {
        "children": [
            {"id": "n1", "x": 0, "y": 0, "width": 100, "height": 50},
            {"id": "n2", "x": 200, "y": 100, "width": 100, "height": 50},
        ],
        "edges": [
            {
                "id": "e1",
                "sections": [
                    {
                        "startPoint": {"x": 100, "y": 25},
                        "endPoint": {"x": 200, "y": 125},
                        "bendPoints": [{"x": 150, "y": 75}],
                    }
                ],
            }
        ],
    }
    result = render_layout_svg(layout_result, nodes)

    # Verify edge path
    edge_path = result.find(".//svg:path", namespaces={"svg": SVG_NS})
    assert edge_path is not None
    assert edge_path.get("id") == "e1"
    assert edge_path.get("d") == "M 100,25 L 150,75 L 200,125"

    # Verify SVG structure
    assert result.tag == "{%s}svg" % SVG_NS
    assert result.nsmap[None] == SVG_NS
    assert "viewBox" in result.attrib

    # Verify groups
    groups = result.findall("svg:g", namespaces={"svg": SVG_NS})
    assert len(groups) == 2  # edges and nodes groups

    # Verify nodes
    nodes = result.findall('.//svg:g[@id="nodes"]/*', namespaces={"svg": SVG_NS})
    assert len(nodes) == 2

    # Verify node content
    rect = result.find('.//svg:rect[@fill="blue"]', namespaces={"svg": SVG_NS})
    circle = result.find('.//svg:circle[@fill="red"]', namespaces={"svg": SVG_NS})
    assert rect is not None
    assert circle is not None


def test_layout_with_svgs():
    """Test complete layout workflow."""
    # Create test SVGs
    nodes = {
        "n1": etree.fromstring('<svg width="100" height="50"></svg>', parser=parser),
        "n2": etree.fromstring('<svg width="120" height="60"></svg>', parser=parser),
    }
    edges = {"e1": ("n1", "n2")}

    # Test basic layout
    result = layout_with_svgs(nodes, edges)
    assert result.tag == "{%s}svg" % SVG_NS
    assert result.nsmap[None] == SVG_NS

    # Test with options
    options = {"elk.direction": "DOWN"}
    result = layout_with_svgs(nodes, edges, options)
    assert result.tag == "{%s}svg" % SVG_NS
