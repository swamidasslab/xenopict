"""Tests for the declarative API."""

import pytest
from pathlib import Path
import numpy as np
from xenopict.declarative import create_molecules
from xenopict import Xenopict


def test_basic_molecule():
    """Test creating a single molecule."""
    spec = {
        "molecules": {
            "smiles": "CCO"
        }
    }
    xenopicts = create_molecules(spec)
    assert len(xenopicts) == 1
    assert isinstance(xenopicts[0], Xenopict)


def test_invalid_smiles():
    """Test that invalid SMILES raises error."""
    spec = {
        "molecules": {
            "smiles": "invalid"
        }
    }
    with pytest.raises(ValueError, match="Invalid SMILES"):
        create_molecules(spec)


def test_missing_smiles():
    """Test that missing SMILES raises error."""
    spec = {
        "molecules": {}
    }
    with pytest.raises(ValueError, match="Field required"):
        create_molecules(spec)


def test_multiple_molecules():
    """Test creating multiple molecules."""
    spec = {
        "molecules": [
            {
                "smiles": "CCO"
            },
            {
                "smiles": "CCCO"
            }
        ]
    }
    xenopicts = create_molecules(spec)
    assert len(xenopicts) == 2
    assert all(isinstance(x, Xenopict) for x in xenopicts)


def get_atom_coords(mol):
    """Get atom coordinates from a molecule."""
    conf = mol.GetConformer()
    return np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])


def test_alignment_enabled():
    """Test that alignment is enabled by default."""
    spec = {
        "molecules": [
            {
                "smiles": "CCO"  # ethanol
            },
            {
                "smiles": "CCCO"  # propanol
            }
        ]
    }
    xenopicts = create_molecules(spec)
    assert len(xenopicts) == 2
    
    # Get coordinates for both molecules
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)

    # The coordinates of the first molecule should match the last three atoms
    # of the second molecule (the CCO portion)
    assert np.allclose(coords1, coords2[1:], atol=0.1)


def test_alignment_disabled():
    """Test that alignment can be disabled."""
    spec = {
        "molecules": [
            {
                "smiles": "CCO"  # ethanol
            },
            {
                "smiles": "CCCO"  # propanol
            }
        ],
        "align": False
    }
    xenopicts = create_molecules(spec)
    assert len(xenopicts) == 2
    
    # Get coordinates for both molecules
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)
    
    # The coordinates should be different since alignment is disabled
    # Note: This could theoretically fail if RDKit happens to generate
    # the same coordinates by chance, but it's extremely unlikely
    assert not np.allclose(coords1, coords2[1:], atol=0.1)


def test_explicit_alignment_enabled():
    """Test that alignment can be explicitly enabled."""
    spec = {
        "molecules": [
            {
                "smiles": "CCO"  # ethanol
            },
            {
                "smiles": "CCCO"  # propanol
            }
        ],
        "align": True
    }
    xenopicts = create_molecules(spec)
    assert len(xenopicts) == 2
    
    # Get coordinates for both molecules
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)
    # The coordinates of the first molecule should match the last three atoms
    # of the second molecule (the CCO portion)
    assert np.allclose(coords1, coords2[1:], atol=0.1)


def test_from_json_string():
    """Test creating molecules from a JSON string."""
    json_str = '''
    {
        "molecules": [
            {
                "smiles": "CCO"
            },
            {
                "smiles": "CCCO"
            }
        ]
    }
    '''
    xenopicts = create_molecules(json_str)
    assert len(xenopicts) == 2
    assert all(isinstance(x, Xenopict) for x in xenopicts)
    
    # Verify alignment works with JSON input
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)
    assert np.allclose(coords1, coords2[1:], atol=0.1)


def test_invalid_json_string():
    """Test that invalid JSON string raises error."""
    json_str = '''
    {
        "molecules": [
            {
                "smiles": "CCO"
            },
    '''
    with pytest.raises(ValueError, match="Invalid JSON string"):
        create_molecules(json_str)


def test_from_file(tmp_path):
    """Test creating molecules from a JSON file."""
    json_str = '''
    {
        "molecules": [
            {
                "smiles": "CCO"
            },
            {
                "smiles": "CCCO"
            }
        ]
    }
    '''
    path = tmp_path / "molecules.json"
    path.write_text(json_str)
    xenopicts = create_molecules(path)
    assert len(xenopicts) == 2
    assert all(isinstance(x, Xenopict) for x in xenopicts)
    
    # Verify alignment works with file input
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)
    assert np.allclose(coords1, coords2[1:], atol=0.1)


def test_file_not_found():
    """Test that non-existent file raises error."""
    with pytest.raises(ValueError, match="File not found"):
        create_molecules(Path("nonexistent.json"))


def test_invalid_json_file(tmp_path):
    """Test that invalid JSON file raises error."""
    json_str = '''
    {
        "molecules": [
            {
                "smiles": "CCO"
            },
    '''
    path = tmp_path / "invalid.json"
    path.write_text(json_str)
    with pytest.raises(ValueError, match="Invalid JSON in file"):
        create_molecules(path) 