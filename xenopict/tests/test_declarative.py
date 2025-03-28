"""Tests for the declarative API."""

import pytest
from pathlib import Path
from xenopict.declarative import parse
from xenopict.declarative.types import XenopictSpec, MoleculeSpec, MarkSpec
import json
from typing import cast
import numpy as np
from numpy.testing import assert_allclose
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
    spec = XenopictSpec(
        molecules=MoleculeSpec(smiles="CCO", id="mol1", mark=None),
        align=True
    )
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
            MoleculeSpec(smiles="CCO", id="mol1", mark=None),
            MoleculeSpec(smiles="CCCO", id="mol2", mark=None)
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


def test_basic_alignment():
    """Test that molecules are aligned by default."""
    spec = {
        "molecules": [
            {"smiles": "CCO", "id": "mol1"},
            {"smiles": "CCCO", "id": "mol2"}
        ]
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 2

    # Get coordinates for both molecules
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)

    # The coordinates of the first molecule should match the last three atoms
    # of the second molecule (the CCO portion)
    assert_allclose(coords1, coords2[1:], atol=0.1)


def test_alignment_off():
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
    
    # Get coordinates for both molecules
    coords1 = get_atom_coords(xenopicts[0].mol)
    coords2 = get_atom_coords(xenopicts[1].mol)
    
    # The coordinates should be different since alignment is disabled
    # Note: This could theoretically fail if RDKit happens to generate
    # the same coordinates by chance, but it's extremely unlikely
    with pytest.raises(AssertionError):
        assert_allclose(coords1, coords2[1:], atol=0.1)


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
    assert_allclose(coords1, coords2[1:], atol=0.1)


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
    with pytest.raises(AssertionError):
        assert_allclose(coords1, coords2[1:], atol=0.1)


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
    assert_allclose(coords1, coords2[1:], atol=0.1)


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
    assert_allclose(coords1, coords2[1:], atol=0.1)


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
    assert_allclose(coords1, coords2[1:], atol=0.1)


def test_mark_spec_atoms():
    """Test valid atom marking."""
    mark = MarkSpec(
        atoms=[0, 1, 2],
        substructure_atoms=None,
        substructure_bonds=None
    )
    assert mark.atoms == [0, 1, 2]
    assert mark.substructure_atoms is None
    assert mark.substructure_bonds is None


def test_mark_spec_substructure():
    """Test valid substructure marking."""
    mark = MarkSpec(
        atoms=None,
        substructure_atoms=[0, 1, 2],
        substructure_bonds=None
    )
    assert mark.substructure_atoms == [0, 1, 2]
    assert mark.atoms is None
    assert mark.substructure_bonds is None


def test_mark_spec_substructure_with_bonds():
    """Test valid substructure marking with specific bonds."""
    mark = MarkSpec(
        atoms=None,
        substructure_atoms=[0, 1, 2],
        substructure_bonds=[(0, 1), (1, 2)]
    )
    assert mark.substructure_atoms == [0, 1, 2]
    assert mark.substructure_bonds == [(0, 1), (1, 2)]
    assert mark.atoms is None


def test_mark_spec_empty():
    """Test that empty spec is invalid."""
    with pytest.raises(ValueError, match="Must specify at least one marking method"):
        MarkSpec(atoms=None, substructure_atoms=None, substructure_bonds=None)


def test_mark_spec_mixing_methods():
    """Test that mixing marking methods is invalid."""
    with pytest.raises(ValueError, match="Cannot mix marking methods"):
        MarkSpec(
            atoms=[0],
            substructure_atoms=[1, 2],
            substructure_bonds=None
        )


def test_mark_spec_bonds_without_atoms():
    """Test that specifying bonds without atoms is invalid."""
    with pytest.raises(ValueError, match="Cannot specify 'substructure_bonds' without 'substructure_atoms'"):
        MarkSpec(
            atoms=None,
            substructure_atoms=None,
            substructure_bonds=[(0, 1)]
        )


def test_mark_spec_negative_atoms():
    """Test that negative atom indices are invalid."""
    with pytest.raises(ValueError, match="Atom indices must be non-negative"):
        MarkSpec(
            atoms=[-1, 0, 1],
            substructure_atoms=None,
            substructure_bonds=None
        )


def test_mark_spec_negative_substructure_atoms():
    """Test that negative substructure atom indices are invalid."""
    with pytest.raises(ValueError, match="Substructure atom indices must be non-negative"):
        MarkSpec(
            atoms=None,
            substructure_atoms=[-1, 0, 1],
            substructure_bonds=None
        )


def test_mark_spec_negative_bond_atoms():
    """Test that negative bond indices are invalid."""
    with pytest.raises(ValueError, match="Substructure bond indices must be non-negative"):
        MarkSpec(
            atoms=None,
            substructure_atoms=[0, 1, 2],
            substructure_bonds=[(-1, 1)]
        )


def test_mark_spec_duplicate_bonds():
    """Test that duplicate bonds are invalid."""
    with pytest.raises(ValueError, match="Duplicate bonds specified"):
        MarkSpec(
            atoms=None,
            substructure_atoms=[0, 1, 2],
            substructure_bonds=[(0, 1), (0, 1)]
        )


def test_mark_spec_self_bonds():
    """Test that self-bonds are invalid."""
    with pytest.raises(ValueError, match="Bond cannot connect an atom to itself"):
        MarkSpec(
            atoms=None,
            substructure_atoms=[0, 1, 2],
            substructure_bonds=[(0, 0)]
        )


def test_mark_spec_from_dict():
    """Test creating MarkSpec from dictionary."""
    data = {
        "atoms": [0, 1, 2]
    }
    mark = MarkSpec.model_validate(data)
    assert mark.atoms == [0, 1, 2]


def test_mark_spec_in_molecule():
    """Test using MarkSpec within MoleculeSpec."""
    mol = MoleculeSpec(
        smiles="CCO",
        id="test-mol",
        mark=MarkSpec(
            atoms=[0, 1],
            substructure_atoms=None,
            substructure_bonds=None
        )
    )
    assert mol.mark is not None
    assert mol.mark.atoms == [0, 1]


def test_mark_spec_in_json():
    """Test parsing MarkSpec from JSON."""
    data = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [0, 1],
                "substructure_bonds": [[0, 1]]
            }
        }
    }
    spec = XenopictSpec.model_validate(data)
    assert isinstance(spec.molecules, MoleculeSpec)
    assert isinstance(spec.molecules.mark, MarkSpec)
    assert spec.molecules.mark.substructure_atoms == [0, 1]
    assert spec.molecules.mark.substructure_bonds == [(0, 1)] 