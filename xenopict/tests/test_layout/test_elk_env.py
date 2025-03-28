"""Test the JavaScript environment for ELK."""

from xenopict.layout.elk import check_js_env, check_elk_loaded

def test_js_env_basic():
    """Test that the JavaScript environment function executes."""
    env = check_js_env()
    assert isinstance(env, dict)

def test_elk_loaded_basic():
    """Test that ELK loaded function executes."""
    loaded = check_elk_loaded()
    assert isinstance(loaded, bool)

def test_js_environment():
    """Test that the JavaScript environment is properly set up."""
    env = check_js_env()
    # Check that we can access the properties
    assert env['hasWindow'] is True
    assert env['hasDocument'] is True
    assert env['hasProcess'] is True
    assert env['hasXHR'] is True

def test_elk_initialization():
    """Test that ELK is properly loaded."""
    assert check_elk_loaded() is True 