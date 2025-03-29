"""
ELK layout engine integration using mini-racer for JavaScript execution.

This module provides a Python interface to the ELK graph layout algorithm
through the elkjs JavaScript library.
"""

import importlib.resources
import json
from dataclasses import dataclass
from enum import Enum
from typing import Any, Dict, List, Optional, TypeVar, cast

import py_mini_racer
from py_mini_racer._objects import JSPromise

T = TypeVar("T")


class LayoutAlgorithm(str, Enum):
    """Available layout algorithms in ELK."""

    LAYERED = "layered"  # Layer-based layout (Sugiyama)
    STRESS = "stress"  # Stress-based layout
    MRTREE = "mrtree"  # Tree layout
    FORCE = "force"  # Force-based layout
    BOX = "box"  # Simple box layout


class LayoutDirection(str, Enum):
    """Available layout directions in ELK."""

    RIGHT = "RIGHT"
    LEFT = "LEFT"
    DOWN = "DOWN"
    UP = "UP"


class EdgeRouting(str, Enum):
    """Available edge routing styles in ELK."""

    ORTHOGONAL = "ORTHOGONAL"
    POLYLINE = "POLYLINE"
    SPLINES = "SPLINES"


@dataclass
class NodeSizeConfig:
    """Configuration for node sizes in the layout."""

    default_width: float = 100
    default_height: float = 100
    min_width: float = 50
    min_height: float = 50


@dataclass
class LayoutSpacing:
    """Spacing configuration for the layout."""

    node_node: float = 20
    edge_edge: float = 10
    edge_node: float = 10
    component_component: float = 20
    edge_edge_between_layers: float = 10
    edge_node_between_layers: float = 10


@dataclass
class LayoutPadding:
    """Padding configuration for the layout."""

    top: float = 25
    left: float = 25
    bottom: float = 25
    right: float = 25

    def __str__(self) -> str:
        """Convert padding to ELK format."""
        return f"top={self.top},left={self.left},bottom={self.bottom},right={self.right}"


@dataclass
class LayoutConfig:
    """Configuration for ELK layout."""

    algorithm: LayoutAlgorithm = LayoutAlgorithm.LAYERED
    direction: LayoutDirection = LayoutDirection.RIGHT
    edge_routing: EdgeRouting = EdgeRouting.ORTHOGONAL
    aspect_ratio: float = 1.6
    node_size: NodeSizeConfig = NodeSizeConfig()
    spacing: LayoutSpacing = LayoutSpacing()
    padding: LayoutPadding = LayoutPadding()


def _build_layout_options(config: LayoutConfig) -> Dict[str, str]:
    """Convert LayoutConfig to ELK options dictionary."""
    return {
        "elk.algorithm": config.algorithm.value,
        "elk.direction": config.direction.value,
        "elk.edgeRouting": config.edge_routing.value,
        "elk.aspectRatio": str(config.aspect_ratio),
        "elk.padding": str(config.padding),
        "elk.spacing.nodeNode": str(config.spacing.node_node),
        "elk.spacing.edgeEdge": str(config.spacing.edge_edge),
        "elk.spacing.edgeNode": str(config.spacing.edge_node),
        "elk.spacing.componentComponent": str(config.spacing.component_component),
        "elk.layered.spacing.edgeEdgeBetweenLayers": str(config.spacing.edge_edge_between_layers),
        "elk.layered.spacing.edgeNodeBetweenLayers": str(config.spacing.edge_node_between_layers),
    }


def _apply_node_sizes(graph: Dict[str, Any], node_size: NodeSizeConfig) -> None:
    """Apply default and minimum node sizes to the graph."""
    if "children" in graph:
        for node in graph["children"]:
            if "width" not in node:
                node["width"] = node_size.default_width
            if "height" not in node:
                node["height"] = node_size.default_height
            # Ensure minimum size
            node["width"] = max(node["width"], node_size.min_width)
            node["height"] = max(node["height"], node_size.min_height)


def elk_to_layout(
    graph: Dict[str, Any],
    config: Optional[LayoutConfig] = None,
    extra_options: Optional[Dict[str, str]] = None,
) -> Dict[str, Any]:
    """
    Apply ELK layout to a graph.

    Args:
        graph: A dictionary representing the graph in ELK JSON format
        config: Optional layout configuration. If not provided, default values will be used.
        extra_options: Additional ELK options to override or extend the configuration.

    Returns:
        A dictionary containing the laid out graph with position information

    Example:
        >>> graph = {
        ...     "id": "root",
        ...     "children": [
        ...         {"id": "n1"},  # Width and height will use defaults
        ...         {"id": "n2", "width": 30, "height": 30}  # Explicit size
        ...     ],
        ...     "edges": [
        ...         {"id": "e1", "sources": ["n1"], "targets": ["n2"]}
        ...     ]
        ... }
        >>> config = LayoutConfig(
        ...     algorithm=LayoutAlgorithm.LAYERED,
        ...     direction=LayoutDirection.RIGHT,
        ...     node_size=NodeSizeConfig(default_width=150, default_height=100)
        ... )
        >>> result = elk_to_layout(graph, config)
        >>> isinstance(result["children"][0]["x"], (int, float))
        True
    """
    # Use default config if none provided
    layout_config = config or LayoutConfig()

    # Build options dictionary from config
    options = _build_layout_options(layout_config)

    # Add any extra options
    if extra_options:
        options.update(extra_options)

    # Pre-process the graph to set node sizes
    _apply_node_sizes(graph, layout_config.node_size)

    # Convert input to JSON
    graph_json = json.dumps(graph)
    options_json = json.dumps(options)

    # Call the global layout function and convert result to JSON string
    r = _ctx.eval(f"""
    (async () => {{
        const graph = {graph_json};
        const options = {options_json};
        const result = await layout(graph, options);
        return JSON.stringify(result);
    }})()
    """)

    if isinstance(r, JSPromise):
        r = r.get(timeout=5)  # Wait for the Promise to resolve with a 5-second timeout
    return json.loads(str(r))


def _eval_async_js(js_code: str) -> Any:
    """
    Evaluate async JavaScript code and return the result.
    Handles promise resolution in both async and sync contexts.
    """
    promise = cast(JSPromise, _ctx.eval(js_code))
    return promise.get()


# Initialize V8 context with ELK
_ctx = py_mini_racer.MiniRacer()

# Set up browser globals
_ctx.eval("""
var globalThis = globalThis || this;
var window = globalThis;
var document = { currentScript: { src: '' } };
var process = { env: {} };
Error.stackTraceLimit = 64;

// Initialize required ELK globals
var XMLHttpRequest = function() {
    this.open = function() {};
    this.send = function() {};
    this.setRequestHeader = function() {};
};

// Test function to verify JavaScript environment
function testEnv() {
    return {
        hasWindow: typeof window !== 'undefined',
        hasDocument: typeof document !== 'undefined',
        hasProcess: typeof process !== 'undefined',
        hasXHR: typeof XMLHttpRequest !== 'undefined'
    };
}

// Override Atomics.waitAsync to use Promise.resolve
Atomics.waitAsync = function() {
    return {
        value: new Promise((resolve) => {
            resolve();
        })
    };
};
""")

# Load ELK library
_elk_js = importlib.resources.files("xenopict.layout").joinpath("js/elk.js").read_text()
_ctx.eval(_elk_js)

# Initialize ELK and expose it globally
_ctx.eval("""
// Create ELK instance with default options
const defaultOptions = {
    'elk.algorithm': 'layered',
    'elk.spacing.nodeNode': '20',
    'elk.layered.spacing.nodeNodeBetweenLayers': '20'
};

// Create a simple layout function
async function layout(graph, options = {}) {
    // Merge options
    const layoutOptions = { ...defaultOptions, ...options };
    
    try {
        // Apply layout using ELK's default export
        return await ELK.layout(graph, layoutOptions);
    } catch (error) {
        console.error('ELK layout error:', error);
        throw error;
    }
}

// Helper function to escape XML strings
function escapeXml(str) {
    return str.replace(/[<>&'"]/g, function(c) {
        switch (c) {
            case '<': return '&lt;';
            case '>': return '&gt;';
            case '&': return '&amp;';
            case "'": return '&apos;';
            case '"': return '&quot;';
        }
    });
}

// Create a simple SVG rendering function
function toSvg(graph, width, height) {
    let svg = '<?xml version="1.0" encoding="UTF-8"?>\n';
    svg += '<svg width="' + width + '" height="' + height + '" ';
    svg += 'xmlns="http://www.w3.org/2000/svg">\n';
    svg += '    <defs>\n';
    svg += '        <marker id="arrowhead" ';
    svg += 'markerWidth="10" markerHeight="7" ';
    svg += 'refX="9" refY="3.5" orient="auto">\n';
    svg += '            <polygon points="0 0, 10 3.5, 0 7" fill="#000"/>\n';
    svg += '        </marker>\n';
    svg += '    </defs>\n';

    for (let node of graph.children || []) {
        // Draw node rectangle
        svg += '    <rect ';
        svg += 'x="' + node.x + '" ';
        svg += 'y="' + node.y + '" ';
        svg += 'width="' + node.width + '" ';
        svg += 'height="' + node.height + '" ';
        svg += 'fill="none" stroke="#000" stroke-width="1"/>\n';

        // Draw node text
        svg += '    <text ';
        svg += 'x="' + (node.x + node.width/2) + '" ';
        svg += 'y="' + (node.y + node.height/2) + '" ';
        svg += 'text-anchor="middle" dominant-baseline="middle">';
        svg += escapeXml(node.id) + '</text>\n';
    }

    // Draw edges
    for (let edge of graph.edges || []) {
        let path = '';
        if (edge.sections && edge.sections.length > 0) {
            const section = edge.sections[0];
            path += 'M ' + section.startPoint.x + ' ' + section.startPoint.y + ' ';
            
            if (section.bendPoints) {
                for (let point of section.bendPoints) {
                    path += 'L ' + point.x + ' ' + point.y + ' ';
                }
            }
            
            path += 'L ' + section.endPoint.x + ' ' + section.endPoint.y;
        }
        
        if (path) {
            svg += '    <path ';
            svg += 'd="' + path + '" ';
            svg += 'fill="none" stroke="#000" ';
            svg += 'stroke-width="1" marker-end="url(#arrowhead)"/>\n';
        }
    }
    
    svg += "</svg>";
    return svg;
}

// Create functions to get layout options and algorithms
function getLayoutOptions() {
    return {
        'elk.algorithm': 'layered',
        'elk.spacing.nodeNode': '20',
        'elk.layered.spacing.nodeNodeBetweenLayers': '20',
        'elk.direction': 'RIGHT',
        'elk.aspectRatio': '1.6',
        'elk.padding': 'top=25,left=25,bottom=25,right=25',
        'elk.edgeRouting': 'ORTHOGONAL',
        'elk.layered.crossingMinimization.strategy': 'LAYER_SWEEP',
        'elk.layered.nodePlacement.strategy': 'BRANDES_KOEPF',
        'elk.layered.spacing.edgeEdgeBetweenLayers': '10',
        'elk.layered.spacing.edgeNodeBetweenLayers': '10',
        'elk.spacing.componentComponent': '20',
        'elk.spacing.edgeEdge': '10',
        'elk.spacing.edgeNode': '10'
    };
}

function getLayoutAlgorithms() {
    return [
        {
            'id': 'layered',
            'name': 'Layer-Based',
            'description': 'Layer-based layout algorithm (Sugiyama)'
        },
        {
            'id': 'stress',
            'name': 'Stress-Based',
            'description': 'Stress-based layout algorithm'
        },
        {
            'id': 'mrtree',
            'name': 'Mr. Tree',
            'description': 'Tree layout algorithm'
        },
        {
            'id': 'force',
            'name': 'Force',
            'description': 'Force-based layout algorithm'
        },
        {
            'id': 'box',
            'name': 'Box Layout',
            'description': 'Simple box layout algorithm'
        }
    ];
}

// Expose functions globally
globalThis.layout = layout;
globalThis.toSvg = toSvg;
globalThis.getLayoutOptions = getLayoutOptions;
globalThis.getLayoutAlgorithms = getLayoutAlgorithms;
""")


def check_js_env() -> Dict[str, bool]:
    """Test the JavaScript environment setup."""
    result = _ctx.eval("""
    (() => {
        const env = testEnv();
        return JSON.stringify(env);
    })();
    """)
    return json.loads(str(result))


def check_elk_loaded() -> bool:
    """Test if ELK is properly loaded."""
    result = _ctx.eval("""
    (() => {
        const loaded = typeof layout !== 'undefined';
        return JSON.stringify(loaded);
    })();
    """)
    return json.loads(str(result))


def elk_to_svg(
    graph: Dict[str, Any],
    config: Optional[LayoutConfig] = None,
    extra_options: Optional[Dict[str, str]] = None,
) -> str:
    """
    Convert an ELK graph directly to SVG, applying layout in the process.

    Args:
        graph: A dictionary representing the graph in ELK JSON format
        config: Optional layout configuration. If not provided, default values will be used.
        extra_options: Additional ELK options to override or extend the configuration.

    Returns:
        A string containing the SVG representation of the laid out graph

    Example:
        >>> graph = {
        ...     "id": "root",
        ...     "children": [
        ...         {"id": "n1", "width": 30, "height": 30},
        ...         {"id": "n2", "width": 30, "height": 30}
        ...     ],
        ...     "edges": [
        ...        {"id": "e1", "sources": ["n1"], "targets": ["n2"]}
        ...     ]
        ... }
        >>> config = LayoutConfig(
        ...     algorithm=LayoutAlgorithm.LAYERED,
        ...     direction=LayoutDirection.RIGHT
        ... )
        >>> svg = elk_to_svg(graph, config)
        >>> svg.startswith('<?xml') or svg.startswith('<svg')
        True
    """
    # First apply the layout, then convert to SVG
    laid_out_graph = elk_to_layout(graph, config, extra_options)
    return layout_to_svg(laid_out_graph)


def layout_to_svg(layout: Dict[str, Any]) -> str:
    """
    Convert a laid out graph to SVG without applying layout again.

    Args:
        layout: A dictionary representing an already laid out graph with position information

    Returns:
        A string containing the SVG representation of the laid out graph

    Example:
        >>> layout = {
        ...     "id": "root",
        ...     "children": [
        ...         {"id": "n1", "width": 30, "height": 30, "x": 0, "y": 0},
        ...         {"id": "n2", "width": 30, "height": 30, "x": 100, "y": 0}
        ...     ],
        ...     "edges": [
        ...        {"id": "e1", "sources": ["n1"], "targets": ["n2"]}
        ...     ]
        ... }
        >>> svg = layout_to_svg(layout)
        >>> svg.startswith('<?xml') or svg.startswith('<svg')
        True
    """
    # Convert to SVG using the renderer
    r = _ctx.eval(f"""
    (() => {{
        const graph = {json.dumps(layout)};
        return toSvg(graph);
    }})()
    """)

    return str(r)


def get_layout_options() -> Dict[str, str]:
    """
    Get available ELK layout options.

    Returns:
        A dictionary containing available layout options and their default values.
        Use LayoutConfig to configure these options in a type-safe way.

    Example:
        >>> options = get_layout_options()
        >>> isinstance(options, dict)
        True
        >>> 'elk.algorithm' in options
        True
    """
    r = _ctx.eval("""
    (() => {
        return JSON.stringify(getLayoutOptions());
    })()
    """)
    return json.loads(str(r))


def get_layout_algorithms() -> List[Dict[str, str]]:
    """
    Get available ELK layout algorithms.

    Returns:
        A list of dictionaries containing layout algorithms and their metadata.
        Use LayoutAlgorithm enum for type-safe algorithm selection.

    Example:
        >>> algorithms = get_layout_algorithms()
        >>> isinstance(algorithms, list)
        True
        >>> all(isinstance(a, dict) for a in algorithms)
        True
        >>> all('id' in a for a in algorithms)
        True
    """
    r = _ctx.eval("""
    (() => {
        return JSON.stringify(getLayoutAlgorithms());
    })()
    """)
    return json.loads(str(r))
