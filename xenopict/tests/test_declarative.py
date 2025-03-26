"""Tests for xenopict's declarative API."""

import pytest
from xenopict.declarative import from_spec
from xenopict import Xenopict

def test_basic_smiles():
    """Test basic SMILES molecule specification."""
    spec = {
        "molecule": {
            "smiles": "CC(=O)O"  # Acetic acid
        }
    }
    result = from_spec(spec)
    assert len(result) == 1
    assert isinstance(result[0], Xenopict)
    assert str(result[0]) != ""  # Should produce non-empty SVG

def test_basic_smarts():
    """Test basic SMARTS pattern specification."""
    spec = {
        "molecule": {
            "smarts": "[CH3][CX3](=O)[OX2H1]"  # Carboxylic acid pattern
        }
    }
    result = from_spec(spec)
    assert len(result) == 1
    assert isinstance(result[0], Xenopict)
    assert str(result[0]) != ""  # Should produce non-empty SVG

def test_invalid_smiles():
    """Test that invalid SMILES raises appropriate error."""
    spec = {
        "molecule": {
            "smiles": "this_is_not_valid_smiles"
        }
    }
    with pytest.raises(ValueError, match="Invalid SMILES"):
        from_spec(spec)

def test_invalid_smarts():
    """Test that invalid SMARTS raises appropriate error."""
    spec = {
        "molecule": {
            "smarts": "this_is_not_valid_smarts"
        }
    }
    with pytest.raises(ValueError, match="Invalid SMARTS"):
        from_spec(spec)

def test_mutually_exclusive_smiles_smarts():
    """Test that providing both SMILES and SMARTS raises error."""
    spec = {
        "molecule": {
            "smiles": "CC(=O)O",
            "smarts": "[CH3][CX3](=O)[OX2H1]"
        }
    }
    with pytest.raises(ValueError, match="Only one of 'smiles' or 'smarts' should be provided"):
        from_spec(spec)

def test_missing_smiles_and_smarts():
    """Test that providing neither SMILES nor SMARTS raises error."""
    spec = {
        "molecule": {
            "backbone_color": "#FF0000"
        }
    }
    with pytest.raises(ValueError, match="Either 'smiles' or 'smarts' must be provided"):
        from_spec(spec)

def test_multiple_molecules():
    """Test handling multiple molecules with mix of SMILES and SMARTS."""
    spec = {
        "molecule": [
            {
                "smiles": "CC(=O)O",
                "backbone_color": "#FF0000"
            },
            {
                "smarts": "[CH3][CX3](=O)[OX2H1]",
                "backbone_color": "#0000FF"
            }
        ]
    }
    result = from_spec(spec)
    assert len(result) == 2
    assert all(isinstance(x, Xenopict) for x in result)
    assert all(str(x) != "" for x in result)  # All should produce non-empty SVGs 