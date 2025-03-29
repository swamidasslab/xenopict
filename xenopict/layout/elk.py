"""
ELK layout engine integration using mini-racer for JavaScript execution.

This module provides a Python interface to the ELK graph layout algorithm
through the elkjs JavaScript library.
"""

import importlib.resources
import json
from typing import Any, Dict, List, Optional, TypeVar, cast

import py_mini_racer
from py_mini_racer._objects import JSPromise

T = TypeVar("T")


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
        const loaded = typeof ELK !== 'undefined';
        return JSON.stringify(loaded);
    })();
    """)
    return json.loads(str(result))


def layout(graph: Dict[str, Any], options: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """
    Apply ELK layout to a graph.

    Args:
        graph: A dictionary representing the graph in ELK JSON format
        options: Optional layout options to override defaults

    Returns:
        A dictionary containing the laid out graph with position information

    Example:
        >>> graph = {
        ...     "id": "root",
        ...     "children": [
        ...         {"id": "n1", "width": 30, "height": 30},
        ...         {"id": "n2", "width": 30, "height": 30}
        ...     ],
        ...     "edges": [
        ...         {"id": "e1", "sources": ["n1"], "targets": ["n2"]}
        ...     ]
        ... }
        >>> result = layout(graph)
        >>> isinstance(result["children"][0]["x"], (int, float))
        True
    """
    # Convert input to JSON
    graph_json = json.dumps(graph)
    options_json = json.dumps(options or {})

    # Call the global layout function and convert result to JSON string
    r = _ctx.eval(f"""
    (async () => {{
        const graph = JSON.parse('{graph_json}');
        const options = JSON.parse('{options_json}');
        const result = await layout(graph, options);
        return JSON.stringify(result);
    }})()
    """)

    if isinstance(r, JSPromise):
        r = r.get(timeout=5)  # Wait for the Promise to resolve with a 5-second timeout
    return json.loads(str(r))


def layout_to_svg(graph: Dict[str, Any], options: Optional[Dict[str, Any]] = None) -> str:
    """
    Apply ELK layout to a graph and convert it to SVG.

    Args:
        graph: A dictionary representing the graph in ELK JSON format
        options: Optional layout options to override defaults

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
        >>> svg = layout_to_svg(graph)
        >>> svg.startswith('<?xml') or svg.startswith('<svg')
        True
    """
    # First apply the layout
    laid_out_graph = layout(graph, options)

    # Call the function with our graph
    result = _ctx.call("elkSvgRenderer.toSvg", laid_out_graph)
    result = _ctx.call("JSON.stringify", result)
    result = json.loads(str(result))
    return result


def get_layout_options() -> List[Dict[str, Any]]:
    """
    Get available ELK layout options.

    Returns:
        A list of dictionaries containing layout options and their metadata
    """
    opt = _ctx.call("ELK.knownLayoutOptions")
    opt = _ctx.call("JSON.stringify", opt)
    opt = json.loads(str(opt))
    return opt


def get_layout_algorithms() -> List[Dict[str, Any]]:
    """
    Get available ELK layout algorithms.

    Returns:
        A list of dictionaries containing layout algorithms and their metadata
    """
    alg = _ctx.call("ELK.knownLayoutAlgorithms")
    alg = _ctx.call("JSON.stringify", alg)
    alg = json.loads(str(alg))
    return alg
