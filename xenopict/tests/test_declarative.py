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
import re


def _debug_svg_structure(xenopict: Xenopict) -> str:
    """Debug helper to print SVG structure."""
    svg = xenopict.to_svg()
    # Add newlines for readability
    svg = re.sub(r'><', '>\n<', svg)
    return svg


def _has_mark_elements(xenopict: Xenopict) -> bool:
    """Check if the SVG contains mark elements."""
    groups = xenopict.svgdom.getElementsByTagName("g")
    for group in groups:
        if group.getAttribute("class") == "mark":
            return True
    return False


def _get_marked_atom_indices(xenopict: Xenopict) -> set[int]:
    """Get the indices of atoms that are marked with circles."""
    # For atom marking, look for circles in the mark group
    mark_groups = [g for g in xenopict.svgdom.getElementsByTagName("g") 
                  if g.getAttribute("class") == "mark"]
    if not mark_groups:
        return set()
    
    marked_atoms = set()
    for mark_group in mark_groups:
        # Look for circles with atom-X class
        for circle in mark_group.getElementsByTagName("circle"):
            class_name = circle.getAttribute("class")
            if class_name and class_name.startswith("atom-"):
                try:
                    atom_idx = int(class_name.split("-")[1])
                    marked_atoms.add(atom_idx)
                except (IndexError, ValueError):
                    continue
    return marked_atoms


def _get_marked_bonds(xenopict: Xenopict) -> set[tuple[int, int]]:
    """Get the pairs of atom indices connected by marked bonds."""
    # For substructure marking, look for paths in the mark group
    mark_groups = [g for g in xenopict.svgdom.getElementsByTagName("g") 
                  if g.getAttribute("class") == "mark"]
    if not mark_groups:
        return set()
    
    marked_bonds = set()
    for mark_group in mark_groups:
        # Look for paths with bond-X atom-Y atom-Z classes
        for path in mark_group.getElementsByTagName("path"):
            class_name = path.getAttribute("class")
            if class_name and "bond-" in class_name and "atom-" in class_name:
                try:
                    # Extract atom indices from class names like "bond-0 atom-0 atom-1"
                    atoms = [int(x.split("-")[1]) for x in class_name.split() if x.startswith("atom-")]
                    if len(atoms) == 2:
                        marked_bonds.add(tuple(sorted(atoms)))  # Sort to ensure consistent ordering
                except (IndexError, ValueError):
                    continue
    
    return marked_bonds


def _count_mark_circles(xenopict: Xenopict) -> int:
    """Count the number of circle elements in the mark layer."""
    svg = xenopict.to_svg()
    # Find all circle elements within the mark group
    circles = re.findall(r'<g class="mark">.*?<circle.*?</g>', svg, re.DOTALL)
    return len(circles)


def _count_mark_paths(xenopict: Xenopict) -> int:
    """Count the number of path elements in the mark layer."""
    svg = xenopict.to_svg()
    # Find all path elements within the mark group
    paths = re.findall(r'<g class="mark">.*?<path.*?</g>', svg, re.DOTALL)
    return len(paths)


def _get_mark_path_d_attrs(xenopict: Xenopict) -> list[str]:
    """Get the d attributes of all path elements in the mark layer."""
    svg = xenopict.to_svg()
    # Extract d attributes from path elements within mark group
    paths = re.findall(r'<g class="mark">.*?<path.*?d="([^"]*)".*?</g>', svg, re.DOTALL)
    return paths


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


def test_parse_single_molecule():
    """Test parsing a single molecule."""
    spec = {
        "molecules": {
            "smiles": "CCO"
        }
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 1
    assert isinstance(xenopicts[0], Xenopict)


def test_parse_multiple_molecules():
    """Test parsing multiple molecules."""
    spec = {
        "molecules": [
            {"smiles": "CCO"},
            {"smiles": "CCCO"}
        ]
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 2
    assert all(isinstance(x, Xenopict) for x in xenopicts)


def test_parse_json_string():
    """Test parsing from a JSON string."""
    json_str = '''
    {
        "molecules": {
            "smiles": "CCO"
        }
    }
    '''
    xenopicts = parse(json_str)
    assert len(xenopicts) == 1
    assert isinstance(xenopicts[0], Xenopict)


def test_parse_invalid_smiles():
    """Test parsing invalid SMILES."""
    spec = {
        "molecules": {
            "smiles": "invalid"
        }
    }
    with pytest.raises(ValueError, match="Invalid SMILES"):
        parse(spec)


def test_parse_invalid_json():
    """Test parsing invalid JSON."""
    json_str = '''
    {
        "molecules": {
            "smiles": "CCO"
        }
    '''  # Missing closing brace
    with pytest.raises(ValueError, match="Invalid JSON"):
        parse(json_str)


def test_parse_with_atom_marking():
    """Test parsing with atom marking."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "atoms": [0, 1]
            }
        }
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 1
    
    # Verify that marks are present
    assert _has_mark_elements(xenopicts[0])
    
    # Verify correct atoms are marked
    marked_atoms = _get_marked_atom_indices(xenopicts[0])
    assert marked_atoms == {0, 1}


def test_parse_with_substructure_marking():
    """Test parsing with substructure marking."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [0, 1]
            }
        }
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 1
    
    # Verify that marks are present
    assert _has_mark_elements(xenopicts[0])
    
    # For substructure marking, verify the bond is marked
    marked_bonds = _get_marked_bonds(xenopicts[0])
    assert marked_bonds == {(0, 1)}  # The C-C bond should be marked


def test_parse_with_substructure_and_bonds_marking():
    """Test parsing with substructure and bonds marking."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [0, 1],
                "substructure_bonds": [(0, 1)]
            }
        }
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 1
    
    # Verify that marks are present
    assert _has_mark_elements(xenopicts[0])
    
    # For substructure marking, verify the bond is marked
    marked_bonds = _get_marked_bonds(xenopicts[0])
    assert marked_bonds == {(0, 1)}


def test_parse_multiple_molecules_with_marking():
    """Test parsing multiple molecules with different marking types."""
    spec = {
        "molecules": [
            {
                "smiles": "CCO",
                "mark": {
                    "atoms": [0, 1]
                }
            },
            {
                "smiles": "CCCO",
                "mark": {
                    "substructure_atoms": [0, 1, 2]
                }
            }
        ]
    }
    xenopicts = parse(spec)
    assert len(xenopicts) == 2
    
    # First molecule should have atom marks
    assert _has_mark_elements(xenopicts[0])
    marked_atoms1 = _get_marked_atom_indices(xenopicts[0])
    assert marked_atoms1 == {0, 1}
    
    # Second molecule should have substructure marks
    assert _has_mark_elements(xenopicts[1])
    marked_bonds2 = _get_marked_bonds(xenopicts[1])
    assert marked_bonds2 == {(0, 1), (1, 2)}  # Both C-C bonds should be marked


def test_parse_with_invalid_atom_indices():
    """Test parsing with invalid atom indices."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "atoms": [-1]  # Invalid negative index
            }
        }
    }
    with pytest.raises(ValueError, match="Atom indices must be non-negative"):
        parse(spec)


def test_parse_with_invalid_substructure_indices():
    """Test parsing with invalid substructure indices."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [-1]  # Invalid negative index
            }
        }
    }
    with pytest.raises(ValueError, match="Substructure atom indices must be non-negative"):
        parse(spec)


def test_parse_with_invalid_bond_indices():
    """Test parsing with invalid bond indices."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [0, 1],
                "substructure_bonds": [(-1, 1)]  # Invalid negative index
            }
        }
    }
    with pytest.raises(ValueError, match="Substructure bond indices must be non-negative"):
        parse(spec)


def test_parse_with_duplicate_bonds():
    """Test parsing with duplicate bonds."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [0, 1],
                "substructure_bonds": [(0, 1), (0, 1)]  # Duplicate bond
            }
        }
    }
    with pytest.raises(ValueError, match="Duplicate bonds specified"):
        parse(spec)


def test_parse_with_self_bond():
    """Test parsing with self-referential bond."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [0, 1],
                "substructure_bonds": [(0, 0)]  # Self bond
            }
        }
    }
    with pytest.raises(ValueError, match="Bond cannot connect an atom to itself"):
        parse(spec)


def test_parse_with_bonds_without_atoms():
    """Test parsing with bonds specified but no atoms."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "substructure_bonds": [(0, 1)]  # No atoms specified
            }
        }
    }
    with pytest.raises(ValueError, match="Cannot specify 'substructure_bonds' without 'substructure_atoms'"):
        parse(spec)


def test_parse_with_mixed_marking_methods():
    """Test parsing with both atom and substructure marking."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "atoms": [0],
                "substructure_atoms": [1, 2]
            }
        }
    }
    with pytest.raises(ValueError, match="Cannot mix marking methods"):
        parse(spec)


def test_parse_with_no_marking_method():
    """Test parsing with empty marking specification."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {}
        }
    }
    with pytest.raises(ValueError, match="Must specify at least one marking method"):
        parse(spec)


def test_debug_svg_structure():
    """Test that SVG structure contains expected elements for atom marking."""
    spec = {
        "molecules": {
            "smiles": "CCO",
            "mark": {
                "atoms": [0, 1]
            }
        }
    }
    xenopicts = parse(spec)
    svg = _debug_svg_structure(xenopicts[0])
    
    # Verify basic SVG structure
    assert '<svg' in svg
    assert '</svg>' in svg
    
    # Verify mark group exists within overlay
    assert '<g class="overlay">' in svg
    assert '<g class="mark"' in svg
    
    # Verify circles for marked atoms exist with correct classes
    assert 'class="atom-0"' in svg
    assert 'class="atom-1"' in svg
    
    # Verify basic molecule structure exists
    assert '<g class="lines"' in svg
    assert '<g class="text"' in svg 