"""Tests for xenopict's declarative API."""

import pytest
from xenopict.declarative import from_spec
from xenopict import Xenopict
from rdkit import Chem  # type: ignore

def test_basic_smiles():
    """Test basic SMILES molecule specification."""
    spec = {
        "molecule": {
            "smiles": "CC(=O)O"  # Acetic acid
        }
    }
    result = from_spec(spec)
    assert len(result) == 1
    assert isinstance(result[0], str)  # Should be SVG string
    assert result[0] != ""  # Should produce non-empty SVG

def test_basic_smarts():
    """Test basic SMARTS pattern specification."""
    spec = {
        "molecule": {
            "smarts": "[CH3][CX3](=O)[OX2H1]"  # Carboxylic acid pattern
        }
    }
    result = from_spec(spec)
    assert len(result) == 1
    assert isinstance(result[0], str)  # Should be SVG string
    assert result[0] != ""  # Should produce non-empty SVG

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
    assert all(isinstance(x, str) for x in result)  # All should be SVG strings
    assert all(x != "" for x in result)  # All should produce non-empty SVGs

def test_backbone_color():
    """Test that backbone color is correctly applied in SVG output."""
    test_color = "#FF0000"
    spec = {
        "molecule": {
            "smiles": "CC(=O)O",  # Acetic acid
            "backbone_color": test_color
        }
    }
    result = from_spec(spec)
    assert len(result) == 1
    svg = result[0]
    
    # The color should be applied to the lines group's style attribute
    # The style attribute should contain stroke:#FF0000
    assert 'class="lines"' in svg  # Verify lines group exists
    assert f'stroke:{test_color}' in svg  # Verify color is applied
    
    # Test that color change is visible by comparing with default
    spec_no_color = {
        "molecule": {
            "smiles": "CC(=O)O"  # Same molecule without color
        }
    }
    result_no_color = from_spec(spec_no_color)
    assert len(result_no_color) == 1
    svg_no_color = result_no_color[0]
    
    # The SVGs should be different
    assert svg != svg_no_color
    # The default SVG should not have our test color
    assert f'stroke:{test_color}' not in svg_no_color

def test_mapids_alignment_smiles():
    """Test alignment using SMILES with atom mapping IDs."""
    spec = {
        "molecule": [
            {
                "id": "ethanol",
                "smiles": "CCO",
                "map": "[CH3:1][CH2:2][OH:3]"
            },
            {
                "id": "ethylamine",
                "smiles": "CCN",
                "map": "[CH3:1][CH2:2][NH2:3]",
                "alignment": {
                    "to": "ethanol",
                    "method": "mapids"
                }
            }
        ]
    }
    result = from_spec(spec)
    assert len(result) == 2
    # TODO: Add more specific assertions about alignment once implemented

def test_mapids_alignment_integers():
    """Test alignment using integer list mapping."""
    spec = {
        "molecule": [
            {
                "id": "ethanol",
                "smiles": "CCO",
                "map": [1, 2, 3]
            },
            {
                "id": "ethylamine",
                "smiles": "CCN",
                "alignment": {
                    "to": "ethanol",
                    "method": "mapids",
                    "map": [1, 2, 3]
                }
            }
        ]
    }
    result = from_spec(spec)
    assert len(result) == 2
    # TODO: Add more specific assertions about alignment once implemented

def test_mapids_alignment_partial():
    """Test alignment with partial mapping (unmapped atoms)."""
    spec = {
        "molecule": [
            {
                "id": "ethanol",
                "smiles": "CCO",
                "map": [1, 2, 0]  # Oxygen unmapped
            },
            {
                "id": "propanol",
                "smiles": "CCCO",
                "alignment": {
                    "to": "ethanol",
                    "method": "mapids",
                    "map": [1, 2, 0, 0]  # Last two carbons unmapped
                }
            }
        ]
    }
    result = from_spec(spec)
    assert len(result) == 2
    # TODO: Add more specific assertions about alignment once implemented

def test_mapids_alignment_override():
    """Test overriding default maps in alignment."""
    spec = {
        "molecule": [
            {
                "id": "mol1",
                "smiles": "CCO",
                "map": "[CH3:1][CH2:2][OH:3]"  # Default mapping
            },
            {
                "id": "mol2",
                "smiles": "CCN",
                "map": "[CH3:1][CH2:2][NH2:3]",  # Default mapping
                "alignment": {
                    "to": "mol1",
                    "method": "mapids",
                    "map": "[CH3:2][CH2:1][NH2:3]",  # Override mapping
                    "to_map": "[CH3:2][CH2:1][OH:3]"  # Override target mapping
                }
            }
        ]
    }
    result = from_spec(spec)
    assert len(result) == 2
    # TODO: Add more specific assertions about alignment once implemented

def test_invalid_map_integers():
    """Test validation of integer list maps."""
    spec = {
        "molecule": {
            "id": "mol",
            "smiles": "CCO",
            "map": [1, -1, 2]  # Negative indices not allowed
        }
    }
    with pytest.raises(ValueError, match="Map list must contain non-negative integers"):
        from_spec(spec)

def test_invalid_map_smiles():
    """Test validation of SMILES map format."""
    spec = {
        "molecule": {
            "id": "mol",
            "smiles": "CCO",
            "map": "CCO"  # Missing atom mapping IDs
        }
    }
    with pytest.raises(ValueError, match="Map SMILES must contain atom mapping IDs"):
        from_spec(spec)

def test_circular_alignment():
    """Test detection of circular alignments."""
    spec = {
        "molecule": [
            {
                "id": "mol1",
                "smiles": "CCO",
                "alignment": {
                    "to": "mol3",
                    "method": "auto"
                }
            },
            {
                "id": "mol2",
                "smiles": "CCN",
                "alignment": {
                    "to": "mol1",
                    "method": "auto"
                }
            },
            {
                "id": "mol3",
                "smiles": "CCC",
                "alignment": {
                    "to": "mol2",
                    "method": "auto"
                }
            }
        ]
    }
    with pytest.raises(ValueError, match="Circular alignment detected"):
        from_spec(spec) 