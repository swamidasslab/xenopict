"""SVG integration with ELK layout engine.

This module provides utilities for working with SVG elements in ELK layouts,
including size extraction, layout rendering, and SVG embedding. It uses lxml
for DOM manipulation and integrates with the ELK layout engine.

Key features:
1. SVG size extraction from width/height attributes or viewBox
2. ELK graph creation from SVG nodes
3. Layout rendering with SVG node embedding
4. Edge routing with straight lines and bend points
"""

import copy
from typing import Any, Dict, Optional, Tuple

from lxml import etree

from .elk import layout

# SVG namespace for creating elements
SVG_NS = "http://www.w3.org/2000/svg"
NSMAP = {None: SVG_NS}  # Default namespace


def get_svg_size(svg_element: etree._Element) -> Tuple[float, float]:
    """
    Extract width and height from an SVG element.

    Args:
        svg_element: An lxml Element representing an SVG

    Returns:
        A tuple of (width, height) in pixels

    Example:
        >>> svg = etree.fromstring('<svg width="100" height="50"></svg>')
        >>> get_svg_size(svg)
        (100.0, 50.0)
    """
    # Get width/height from attributes with default of '0'
    width = svg_element.get("width", default="0")
    height = svg_element.get("height", default="0")

    # Strip units if present and convert to float
    width = float(width.replace("px", ""))
    height = float(height.replace("px", ""))

    # If no explicit size, try viewBox
    if width == 0 or height == 0:
        viewbox = svg_element.get("viewBox", default="")
        if viewbox:
            _, _, width, height = map(float, viewbox.split())

    return width, height


def create_elk_graph(
    nodes: Dict[str, etree._Element], edges: Dict[str, Tuple[str, str]]
) -> Dict[str, Any]:
    """
    Create an ELK graph structure from SVG nodes and edges.

    Args:
        nodes: Dictionary mapping node IDs to their SVG elements
        edges: Dictionary mapping edge IDs to (source_id, target_id) tuples

    Returns:
        An ELK graph structure ready for layout

    Example:
        >>> nodes = {'n1': etree.fromstring('<svg width="100" height="50"></svg>')}
        >>> edges = {'e1': ('n1', 'n2')}
        >>> graph = create_elk_graph(nodes, edges)
        >>> isinstance(graph, dict)
        True
    """
    children = []
    for node_id, svg in nodes.items():
        width, height = get_svg_size(svg)
        children.append({"id": node_id, "width": width, "height": height})

    edge_list = []
    for edge_id, (source, target) in edges.items():
        edge_list.append({"id": edge_id, "sources": [source], "targets": [target]})

    return {"id": "root", "children": children, "edges": edge_list}


def render_layout_svg(
    layout_result: Dict[str, Any],
    nodes: Dict[str, etree._Element],
    padding: float = 10.0,
) -> etree._Element:
    """
    Render an ELK layout result as an SVG, embedding the original SVG nodes.

    Args:
        layout_result: The result from ELK layout
        nodes: Dictionary mapping node IDs to their SVG elements
        padding: Extra space around the layout in pixels

    Returns:
        An SVG element containing the rendered layout with embedded SVGs

    Example:
        >>> nodes = {'n1': etree.fromstring('<svg width="100" height="50"></svg>')}
        >>> layout_result = {'children': [{'id': 'n1', 'x': 0, 'y': 0}]}
        >>> svg = render_layout_svg(layout_result, nodes)
        >>> svg.tag == 'svg'
        True
    """
    # Calculate total bounds
    min_x = min(child["x"] for child in layout_result["children"])
    min_y = min(child["y"] for child in layout_result["children"])
    max_x = max(child["x"] + child["width"] for child in layout_result["children"])
    max_y = max(child["y"] + child["height"] for child in layout_result["children"])

    # Create root SVG with padding
    width = max_x - min_x + 2 * padding
    height = max_y - min_y + 2 * padding
    root = etree.Element(
        "svg",
        attrib={
            "width": str(width),
            "height": str(height),
            "viewBox": f"{-padding} {-padding} {width} {height}",
        },
        nsmap=NSMAP,
    )

    # Add edges first (under nodes)
    if "edges" in layout_result:
        edges_group = etree.SubElement(root, "g", attrib={"id": "edges"}, nsmap=NSMAP)
        for edge in layout_result["edges"]:
            if "sections" in edge:
                # Get edge points
                points = []
                for section in edge["sections"]:
                    points.extend(
                        [
                            (section["startPoint"]["x"], section["startPoint"]["y"]),
                            (section["endPoint"]["x"], section["endPoint"]["y"]),
                        ]
                    )
                    if "bendPoints" in section:
                        for point in section["bendPoints"]:
                            points.append((point["x"], point["y"]))

                # Create path
                path_data = f"M {points[0][0]},{points[0][1]}"
                for x, y in points[1:]:
                    path_data += f" L {x},{y}"

                etree.SubElement(
                    edges_group,
                    "path",
                    attrib={
                        "d": path_data,
                        "stroke": "black",
                        "fill": "none",
                        "id": edge["id"],
                    },
                    nsmap=NSMAP,
                )

    # Add nodes with their SVGs
    nodes_group = etree.SubElement(root, "g", attrib={"id": "nodes"}, nsmap=NSMAP)
    for child in layout_result["children"]:
        node_id = child["id"]
        if node_id in nodes:
            # Create group for node
            node_group = etree.SubElement(
                nodes_group,
                "g",
                attrib={"transform": f"translate({child['x']},{child['y']})"},
                nsmap=NSMAP,
            )

            # Clone and embed the original SVG content
            svg = nodes[node_id]
            for elem in svg.getchildren():
                node_group.append(copy.deepcopy(elem))

    return root


def layout_with_svgs(
    nodes: Dict[str, etree._Element],
    edges: Dict[str, Tuple[str, str]],
    options: Optional[Dict[str, Any]] = None,
) -> etree._Element:
    """
    Create a layout from SVG nodes and edges, returning a single SVG.

    This is a convenience function that combines create_elk_graph,
    layout, and render_layout_svg into a single operation.

    Args:
        nodes: Dictionary mapping node IDs to their SVG elements
        edges: Dictionary mapping edge IDs to (source_id, target_id) tuples
        options: Optional ELK layout options

    Returns:
        An SVG element containing the complete layout

    Example:
        >>> nodes = {'n1': etree.fromstring('<svg width="100" height="50"></svg>')}
        >>> edges = {}
        >>> svg = layout_with_svgs(nodes, edges)
        >>> svg.tag == 'svg'
        True
    """
    graph = create_elk_graph(nodes, edges)
    layout_result = layout(graph, options)
    return render_layout_svg(layout_result, nodes)
