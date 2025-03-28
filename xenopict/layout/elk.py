"""
ELK layout engine integration using mini-racer for JavaScript execution.

This module provides a Python interface to the ELK graph layout algorithm
through the elkjs JavaScript library.
"""

import os
import json
import time
import asyncio
from pathlib import Path
from typing import Any, Dict, List, Optional, Union, cast

import py_mini_racer
from py_mini_racer._objects import JSPromise

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
_elk_js_path = Path(__file__).parent / "js" / "elk.js"
with open(_elk_js_path, "r") as f:
    _elk_js = f.read()
    _ctx.eval(_elk_js)

# Initialize ELK
_ctx.eval("const elk = new ELK();")

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

async def layout(graph: Dict[str, Any], options: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
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
        >>> result = asyncio.run(layout(graph))
        >>> isinstance(result["children"][0]["x"], (int, float))
        True
    """
    if options:
        layout_options = json.dumps(options)
        _ctx.eval(f"elk.defaultLayoutOptions = {layout_options};")
    
    # Convert graph to JSON and run layout
    graph_json = json.dumps(graph)
    promise = cast(JSPromise, _ctx.eval(f"""
    (async () => {{
        try {{
            const result = await elk.layout(JSON.parse('{graph_json}'));
            return JSON.stringify(result);
        }} catch (error) {{
            throw new Error('ELK layout failed: ' + error.message);
        }}
    }})();
    """))
    result = await promise
    return json.loads(str(result))

async def get_layout_options() -> List[Dict[str, Any]]:
    """
    Get available ELK layout options.
    
    Returns:
        A list of dictionaries containing all available layout options and their metadata
    """
    promise = cast(JSPromise, _ctx.eval("""
    (async () => {
        try {
            const options = await elk.knownLayoutOptions();
            return JSON.stringify(options);
        } catch (error) {
            throw new Error('Failed to get layout options: ' + error.message);
        }
    })();
    """))
    result = await promise
    return json.loads(str(result))

async def get_layout_algorithms() -> List[Dict[str, Any]]:
    """
    Get available ELK layout algorithms.
    
    Returns:
        A list of dictionaries containing all available layout algorithms and their metadata
    """
    promise = cast(JSPromise, _ctx.eval("""
    (async () => {
        try {
            const algorithms = await elk.knownLayoutAlgorithms();
            return JSON.stringify(algorithms);
        } catch (error) {
            throw new Error('Failed to get layout algorithms: ' + error.message);
        }
    })();
    """))
    result = await promise
    return json.loads(str(result)) 