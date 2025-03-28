"""Tests for basic mini-racer functionality."""

import json

from py_mini_racer import MiniRacer


def test_miniracer_basic():
    """Test basic JavaScript evaluation with mini-racer."""
    ctx = MiniRacer()
    result = ctx.eval("6 * 7")
    assert result == 42


def test_miniracer_json():
    """Test JSON object handling with mini-racer."""
    ctx = MiniRacer()
    result = ctx.eval("""
        var obj = {
            number: 42,
            string: "hello",
            array: [1, 2, 3]
        };
        JSON.stringify(obj);
    """)
    # Parse the JSON string to get a Python dict
    result = json.loads(result)
    assert result == {"number": 42, "string": "hello", "array": [1, 2, 3]}


def test_miniracer_function():
    """Test JavaScript function execution with mini-racer."""
    ctx = MiniRacer()
    ctx.eval("""
        function add(a, b) {
            return a + b;
        }
    """)
    result = ctx.eval("add(40, 2)")
    assert result == 42
