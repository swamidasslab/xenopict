"""Tests for the declarative API."""

import pytest
from pathlib import Path
from xenopict.declarative import parse
from xenopict.declarative.types import XenopictSpec, MoleculeSpec
import json
import numpy as np
from xenopict import Xenopict


def test_single_molecule_dict():
    """Test creating a single molecule from a dict."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "id": "mol1"
        }
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 1
    assert xenopicts[0].mol is not None


def test_single_molecule_spec():
    """Test creating a single molecule from a XenopictSpec."""
    spec = XenopictSpec(molecules=MoleculeSpec(smiles="CCO", id="mol1"), align=True)
    xenopicts = parse(spec)
    assert len(xenopicts) == 1
    assert xenopicts[0].mol is not None


def test_multiple_molecules_dict():
    """Test creating multiple molecules from a dict."""
    spec = {
        "molecules": [
            {"smiles": "CCO", "id": "mol1"},
            {"smiles": "CCCO", "id": "mol2"}
        ]
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 2
    assert all(x.mol is not None for x in xenopicts)


def test_multiple_molecules_spec():
    """Test creating multiple molecules from a XenopictSpec."""
    spec = XenopictSpec(
        molecules=[
            MoleculeSpec(smiles="CCO", id="mol1"),
            MoleculeSpec(smiles="CCCO", id="mol2")
        ],
        align=True
    )
    xenopicts = parse(spec)
    assert len(xenopicts) == 2
    assert all(x.mol is not None for x in xenopicts)


def test_json_string():
    """Test creating molecules from a JSON string."""
    json_str = '''
    {
        "molecules": [
            {
                "smiles": "CCO",
                "id": "mol1"
            },
            {
                "smiles": "CCCO",
                "id": "mol2"
            }
        ]
    }
    '''
    xenopicts = parse(json_str)
    assert len(xenopicts) == 2
    assert all(x.mol is not None for x in xenopicts)


def test_json_file(tmp_path):
    """Test creating molecules from a JSON file."""
    spec = {
        "molecules": [
            {"smiles": "CCO", "id": "mol1"},
            {"smiles": "CCCO", "id": "mol2"}
        ]
    }
    json_file = tmp_path / "molecules.json"
    with open(json_file, "w") as f:
        json.dump(spec, f)
    
    xenopicts = parse(json_file)
    assert len(xenopicts) == 2
    assert all(x.mol is not None for x in xenopicts)


def test_invalid_smiles():
    """Test error handling for invalid SMILES."""
    spec = {
        "molecules": {
            "smiles": "invalid",
            "id": "mol1"
        }
    }
    with pytest.raises(ValueError, match="Invalid SMILES"):
        parse(spec)


def test_invalid_json_string():
    """Test error handling for invalid JSON string."""
    json_str = '''
    {
        "molecules": [
            {
                "smiles": "CCO",
                "id": "mol1"
            },
    '''
    with pytest.raises(ValueError, match="Invalid JSON string"):
        parse(json_str)


def test_missing_file():
    """Test error handling for missing file."""
    with pytest.raises(ValueError, match="File not found"):
        parse(Path("nonexistent.json"))


def test_invalid_json_file(tmp_path):
    """Test error handling for invalid JSON file."""
    json_file = tmp_path / "invalid.json"
    with open(json_file, "w") as f:
        f.write("invalid json")
    
    with pytest.raises(ValueError, match="Invalid JSON in file"):
        parse(json_file)


def test_invalid_input_type():
    """Test error handling for invalid input type."""
    with pytest.raises(TypeError, match="Unsupported input type"):
        parse(42)  # type: ignore


def test_alignment_enabled():
    """Test that molecules are aligned by default."""
    spec = {
        "molecules": [
            {"smiles": "CCO", "id": "mol1"},
            {"smiles": "CCCO", "id": "mol2"}
        ]
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 2
    # TODO: Add more specific alignment checks


def test_alignment_disabled():
    """Test that alignment can be disabled."""
    spec = {
        "molecules": [
            {"smiles": "CCO", "id": "mol1"},
            {"smiles": "CCCO", "id": "mol2"}
        ],
        "align": False
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 2
    # TODO: Add more specific alignment checks


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
    xenopicts = parse(spec)
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
    xenopicts = parse(spec)
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
    xenopicts = parse(spec)
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
    xenopicts = parse(json_str)
    assert len(xenopicts) == 2
    assert all(isinstance(x, Xenopict) for x in xenopicts)
    
    # Verify alignment works with JSON input
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)
    assert np.allclose(coords1, coords2[1:], atol=0.1)


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
    xenopicts = parse(path)
    assert len(xenopicts) == 2
    assert all(isinstance(x, Xenopict) for x in xenopicts)
    
    # Verify alignment works with file input
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)
    assert np.allclose(coords1, coords2[1:], atol=0.1) 