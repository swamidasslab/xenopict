"""Tests for molecular alignment functionality."""

from rdkit import Chem
import pytest
from rdkit.Chem import rdDepictor

from xenopict.alignment import (
    align_to_template,
    align_to_template_manual,
    align_to_template_with_indices,
    align_to_template_by_mapids,
    auto_alignment,
)


def test_align_to_template_no_match():
    """Fails when no substructure match."""
    template = Chem.MolFromSmiles("O")  # oxygen
    mol = Chem.MolFromSmiles("c1ccccc1")  # benzene
    with pytest.raises(ValueError, match="No alignment found"):
        align_to_template(mol, template)


@pytest.mark.parametrize(
    "mol_smiles, template_smiles",
    [
        ("CCCc1ncccc1", "CCc1ncccc1"),
        ("c1cccc(C)c1O", "C1C(C)=CC=CC1=O"),
        ("CC(=O)c1ccccc1", "CC(O)c1ccccc1"),
        ("c1ccccc1CC(=O)O", "c1ccccc1CCO"),
        ("Cc1cccc(O)c1", "CC1=CC=CC(=O)C1"),
        ("c1ccccc1C(=O)O", "c1ccccc1CO"),
    ],
    ids=[
        "2-methylpyridine->2-ethylpyridine",
        "methylphenol->methylcyclohexenone",
        "acetophenone->methylbenzyl alcohol",
        "phenylacetic acid->phenethyl alcohol",
        "methylphenol->methylcyclohexenone",
        "benzoic acid->benzyl alcohol",
    ],
)
def test_align_to_template(template_smiles, mol_smiles):
    """Basic template alignment. This set of tests requires all molecules are uniquely alignable."""

    mol = Chem.MolFromSmiles(mol_smiles)
    template = Chem.MolFromSmiles(template_smiles)

    assert mol
    assert template
    # check that we can find an alignment
    aligned_atoms = auto_alignment(mol, template).aligned_atoms
    assert aligned_atoms, "No alignment found"

    aligned = align_to_template(mol, template)

    assert aligned is mol  # Should modify in place

    # check that the alignment is correct
    _assert_coords_close(mol, template, aligned_atoms)

    # now reverse the alignment as a sanity check
    mol = Chem.MolFromSmiles(mol_smiles)
    template = Chem.MolFromSmiles(template_smiles)

    aligned_atoms = auto_alignment(template, mol).aligned_atoms
    assert aligned_atoms, "No alignment found"
    aligned = align_to_template(template, mol)
    _assert_coords_close(template, mol, aligned_atoms)


@pytest.mark.parametrize(
    "mol_smiles, template_smiles, aligned_atoms",
    [
        (
            "CCc1ccccn1",  # 2-ethylpyridine
            "Cc1ccccn1",  # 2-methylpyridine
            [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5), (6, 6)],
        ),
        (
            "c1ccccc1C(=O)O",
            "c1ccccc1CO",
            [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)],
        ),
    ],
    ids=[
        "2-ethylpyridine->2-methylpyridine",
        "benzoic acid->benzyl alcohol",
    ],
)
def test_align_to_template_manual_pairs(mol_smiles, template_smiles, aligned_atoms):
    """Testing align_to_template_manual with explicit atom pairs"""
    mol = Chem.MolFromSmiles(mol_smiles)
    template = Chem.MolFromSmiles(template_smiles)

    assert mol, f"Failed to parse molecule SMILES: {mol_smiles}"
    assert template, f"Failed to parse template SMILES: {template_smiles}"

    print(f"\nMol atoms: {mol.GetNumAtoms()}")
    print(f"Template atoms: {template.GetNumAtoms()}")
    print(f"Aligned atoms: {aligned_atoms}")

    # Validate atom indices are within bounds
    for mol_idx, template_idx in aligned_atoms:
        assert mol_idx < mol.GetNumAtoms(), (
            f"Molecule atom index {mol_idx} out of range (max {mol.GetNumAtoms() - 1})"
        )
        assert template_idx < template.GetNumAtoms(), (
            f"Template atom index {template_idx} out of range (max {template.GetNumAtoms() - 1})"
        )
        assert mol_idx >= 0, f"Negative molecule atom index {mol_idx}"
        assert template_idx >= 0, f"Negative template atom index {template_idx}"

    aligned = align_to_template_manual(mol, template, aligned_atoms)
    assert aligned is mol  # Should modify in place

    # Check that the alignment is correct
    _assert_coords_close(mol, template, aligned_atoms)

    # Test reverse alignment
    mol = Chem.MolFromSmiles(mol_smiles)
    template = Chem.MolFromSmiles(template_smiles)

    reversed_atoms = [(b, a) for a, b in aligned_atoms]
    aligned = align_to_template_manual(template, mol, reversed_atoms)
    _assert_coords_close(template, mol, reversed_atoms)


def test_align_to_template_manual_invalid_indices():
    """Manual alignment with invalid indices."""
    template = Chem.MolFromSmiles("CCC")  # 3 atoms
    mol = Chem.MolFromSmiles("CC")  # 2 atoms
    rdDepictor.Compute2DCoords(template)
    rdDepictor.Compute2DCoords(mol)

    # Test template index out of bounds
    with pytest.raises(AssertionError, match="Template atom index .* out of range"):
        align_to_template_manual(mol, template, [(0, template.GetNumAtoms())])

    # Test molecule index out of bounds
    with pytest.raises(AssertionError, match="Molecule atom index .* out of range"):
        align_to_template_manual(mol, template, [(mol.GetNumAtoms(), 0)])

    # Test negative indices
    with pytest.raises(AssertionError, match="Negative .* atom index"):
        align_to_template_manual(mol, template, [(-1, 0)])


@pytest.mark.parametrize(
    "mol_smiles,template_smiles,mol_to_template,aligned_atoms",
    [
        ("CC", "CC", [0, 1], [(0, 0), (1, 1)]),  # Simple ethane
        ("CC", "CC", {0: 0, 1: 1}, [(0, 0), (1, 1)]),  # Simple ethane dict
        ("CCC", "CC", [-1, 0, 1], [(1, 0), (2, 1)]),  # Propane to ethane
        ("CCC", "CC", {1: 0, 2: 1}, [(1, 0), (2, 1)]),  # Propane to ethane dict
        ("CC", "CCC", [1, 2], [(0, 1), (1, 2)]),  # Ethane to propane
        ("CC", "CCC", {0: 1, 1: 2}, [(0, 1), (1, 2)]),  # Ethane to propane dict
        (
            "c1ccccc1",
            "c1ccccc1",
            [(i + 1) % 6 for i in range(6)],
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
        ),  # Benzene ring
        (
            "c1ccccc1",
            "c1ccccc1",
            {i: (i + 1) % 6 for i in range(6)},
            [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)],
        ),  # Benzene ring dict
    ],
    ids=[
        "small_to_small_list",
        "small_to_small_dict",
        "large_to_small_list",
        "large_to_small_dict",
        "small_to_large_list",
        "small_to_large_dict",
        "ring_to_ring_list",
        "ring_to_ring_dict",
    ],
)
def test_align_to_template_with_indices_valid(
    mol_smiles, template_smiles, mol_to_template, aligned_atoms
):
    """Test align_to_template_with_indices"""
    mol = Chem.MolFromSmiles(mol_smiles)
    template = Chem.MolFromSmiles(template_smiles)

    assert mol, f"Failed to parse molecule SMILES: {mol_smiles}"
    assert template, f"Failed to parse template SMILES: {template_smiles}"

    aligned = align_to_template_with_indices(mol, template, mol_to_template)
    assert aligned is mol  # Should modify in place

    # Check that the alignment is correct
    _assert_coords_close(mol, template, aligned_atoms)


@pytest.mark.parametrize(
    "mol_smiles,template_smiles,mol_to_template,error_type,error_match",
    [
        (
            "CC",
            "CC",
            [-1, 2],
            ValueError,
            "Reference atom index .* out of range",
        ),  # Invalid template index
        (
            "CC",
            "CC",
            {0: -2},
            ValueError,
            "Reference atom index .* out of range",
        ),  # Invalid template index in dict
        (
            "CC",
            "CC",
            "invalid",
            ValueError,
            "must be a list or dictionary",
        ),  # Invalid type
        (
            "CC",
            "CC",
            [0],
            AssertionError,
            "Length of mol_to_template must equal",
        ),  # Wrong length list
    ],
    ids=["invalid_index", "invalid_dict_index", "invalid_type", "wrong_length"],
)
def test_align_to_template_with_indices_invalid(
    mol_smiles, template_smiles, mol_to_template, error_type, error_match
):
    """Test align_to_template_with_indices with invalid inputs"""
    mol = Chem.MolFromSmiles(mol_smiles)
    template = Chem.MolFromSmiles(template_smiles)

    assert mol, f"Failed to parse molecule SMILES: {mol_smiles}"
    assert template, f"Failed to parse template SMILES: {template_smiles}"

    with pytest.raises(error_type, match=error_match):
        align_to_template_with_indices(mol, template, mol_to_template)


@pytest.mark.parametrize(
    "mol_smiles,template_smiles,mol_map_smiles,template_map_smiles,aligned_atoms",
    [
        (
            "CC",
            "CC",  # Simple 2-atom molecules
            "[CH3:1][CH3:2]",
            "[CH3:1][CH3:2]",
            [(0, 0), (1, 1)],
        ),
        (
            "CCC",
            "CC",  # Larger mol to smaller template
            "[CH3:1][CH2:2][CH3:3]",
            "[CH3:1][CH3:2]",
            [(0, 0), (1, 1)],
        ),
        (
            "CC",
            "CCC",  # Smaller mol to larger template
            "[CH3:1][CH3:2]",
            "[CH3:1][CH2:2][CH3:3]",
            [(0, 0), (1, 1)],
        ),
        (
            "c1ccccc1",
            "c1ccccc1",  # Benzene ring
            "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1",
            "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1",
            [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)],
        ),
    ],
    ids=["small_to_small", "large_to_small", "small_to_large", "ring_to_ring"],
)
def test_align_to_template_by_mapids_valid(
    mol_smiles, template_smiles, mol_map_smiles, template_map_smiles, aligned_atoms
):
    """Test align_to_template_by_mapids with valid inputs"""
    mol = Chem.MolFromSmiles(mol_smiles)
    template = Chem.MolFromSmiles(template_smiles)
    mol_map = Chem.MolFromSmiles(mol_map_smiles)
    template_map = Chem.MolFromSmiles(template_map_smiles)

    assert mol, f"Failed to parse molecule SMILES: {mol_smiles}"
    assert template, f"Failed to parse template SMILES: {template_smiles}"
    assert mol_map, f"Failed to parse mol map SMILES: {mol_map_smiles}"
    assert template_map, f"Failed to parse template map SMILES: {template_map_smiles}"

    aligned = align_to_template_by_mapids(mol, template, mol_map, template_map)
    assert aligned is mol  # Should modify in place

    # Check that the alignment is correct
    _assert_coords_close(mol, template, aligned_atoms)


def GetCoords(mol: Chem.Mol, i: int):
    c = mol.GetConformer(0).GetAtomPosition(i)
    return c.x, c.y, c.z


def _assert_coords_close(
    mol1: Chem.Mol, mol2: Chem.Mol, matched_atoms: list[tuple[int, int]]
):
    for i, j in matched_atoms:
        x1, y1, z1 = GetCoords(mol1, i)
        x2, y2, z2 = GetCoords(mol2, j)
        assert abs(x1 - x2) < 0.01 and abs(y1 - y2) < 0.01 and abs(z1 - z2) < 0.01, (
            f"Atom {i} and atom {j} are not close"
        )
