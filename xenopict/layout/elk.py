"""
ELK layout engine integration using mini-racer for JavaScript execution.

This module provides a Python interface to the ELK graph layout algorithm
through the elkjs JavaScript library.
"""

import json
import sys
from typing import Any, Dict, Optional, Union

import py_mini_racer

if sys.version_info >= (3, 9):
    from importlib import resources
else:
    import importlib_resources as resources

# Initialize V8 context with ELK
_ctx = py_mini_racer.MiniRacer()

# Set up browser globals
_ctx.eval("""
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
""")

# Load ELK library
try:
    # Get reference to the js directory within the package
    if sys.version_info >= (3, 9):
        js_files = resources.files('xenopict.layout.js')
        elk_js = (js_files / 'elk.js').read_text(encoding='utf-8')
    else:
        # Fallback for Python < 3.9
        elk_js = resources.read_text('xenopict.layout.js', 'elk.js')
    
    _ctx.eval(elk_js)
except Exception as e:
    raise RuntimeError(f"Failed to load ELK JavaScript library: {e}")

# Test ELK initialization
_ctx.eval("""
// Simple synchronous test of ELK
const elk = new ELK();
const testGraph = {
    id: "test",
    children: [
        {id: "n1", width: 10, height: 10},
        {id: "n2", width: 10, height: 10}
    ],
    edges: []
};
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
        const loaded = typeof ELK !== 'undefined' && typeof elk !== 'undefined';
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
    if options:
        layout_options = json.dumps(options)
        _ctx.eval(f"elk.defaultLayoutOptions = {layout_options};")
    
    # Convert graph to JSON and run layout
    graph_json = json.dumps(graph)
    result = _ctx.eval(f"""
    (() => {{
        try {{
            const result = elk.layout(JSON.parse('{graph_json}'));
            return JSON.stringify(result);
        }} catch (error) {{
            throw new Error('ELK layout failed: ' + error.message);
        }}
    }})();
    """)
    
    return json.loads(str(result))

def get_layout_options() -> Dict[str, Any]:
    """
    Get available ELK layout options.
    
    Returns:
        A dictionary containing all available layout options and their metadata
    """
    result = _ctx.eval("""
    (() => {
        try {
            const options = elk.knownLayoutOptions();
            return JSON.stringify(options);
        } catch (error) {
            throw new Error('Failed to get layout options: ' + error.message);
        }
    })();
    """)
    return json.loads(str(result))

def get_layout_algorithms() -> Dict[str, Any]:
    """
    Get available ELK layout algorithms.
    
    Returns:
        A dictionary containing all available layout algorithms and their metadata
    """
    result = _ctx.eval("""
    (() => {
        try {
            const algorithms = elk.knownLayoutAlgorithms();
            return JSON.stringify(algorithms);
        } catch (error) {
            throw new Error('Failed to get layout algorithms: ' + error.message);
        }
    })();
    """)
    return json.loads(str(result)) 