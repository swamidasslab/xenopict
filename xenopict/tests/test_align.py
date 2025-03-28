"""Tests for molecular alignment functionality."""

from typing import List, Tuple
from rdkit import Chem
import pytest
from rdkit.Chem import rdDepictor

from xenopict.alignment import (
    align_from_mcs,
    align_from_atom_pairs,
    align_from_indices,
    align_from_mapids,
    auto_align_molecules,
    Alignment,
)


def test_align_from_mcs_no_match():
    """Fails when no substructure match."""
    template_mol = Chem.MolFromSmiles("O")  # oxygen
    source_mol = Chem.MolFromSmiles("c1ccccc1")  # benzene
    with pytest.raises(ValueError, match="No alignment found"):
        align_from_mcs(source_mol, template_mol)


@pytest.mark.parametrize(
    "source_smiles, template_smiles",
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
def test_align_from_mcs(template_smiles, source_smiles):
    """Basic template alignment. This set of tests requires all molecules are uniquely alignable."""

    source_mol = Chem.MolFromSmiles(source_smiles)
    template_mol = Chem.MolFromSmiles(template_smiles)

    assert source_mol
    assert template_mol
    # check that we can find an alignment
    atom_pairs = Alignment.from_mcs(source_mol, template_mol).atom_pairs
    assert atom_pairs, "No alignment found"

    aligned = align_from_mcs(source_mol, template_mol)

    assert aligned is source_mol  # Should modify in place

    # check that the alignment is correct
    _assert_coords_close(source_mol, template_mol, atom_pairs)

    # now reverse the alignment as a sanity check
    source_mol = Chem.MolFromSmiles(source_smiles)
    template_mol = Chem.MolFromSmiles(template_smiles)

    atom_pairs = Alignment.from_mcs(template_mol, source_mol).atom_pairs
    assert atom_pairs, "No alignment found"
    aligned = align_from_mcs(template_mol, source_mol)
    _assert_coords_close(template_mol, source_mol, atom_pairs)


@pytest.mark.parametrize(
    "source_smiles, template_smiles, atom_pairs",
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
def test_align_from_atom_pairs(source_smiles, template_smiles, atom_pairs):
    """Testing align_from_atom_pairs with explicit atom pairs"""
    source_mol = Chem.MolFromSmiles(source_smiles)
    template_mol = Chem.MolFromSmiles(template_smiles)

    assert source_mol, f"Failed to parse source molecule SMILES: {source_smiles}"
    assert template_mol, f"Failed to parse template SMILES: {template_smiles}"

    print(f"\nSource mol atoms: {source_mol.GetNumAtoms()}")
    print(f"Template atoms: {template_mol.GetNumAtoms()}")
    print(f"Atom pairs: {atom_pairs}")

    # Validate atom indices are within bounds
    for source_idx, template_idx in atom_pairs:
        assert source_idx < source_mol.GetNumAtoms(), (
            f"Source molecule atom index {source_idx} out of range (max {source_mol.GetNumAtoms() - 1})"
        )
        assert template_idx < template_mol.GetNumAtoms(), (
            f"Template atom index {template_idx} out of range (max {template_mol.GetNumAtoms() - 1})"
        )
        assert source_idx >= 0, f"Negative source molecule atom index {source_idx}"
        assert template_idx >= 0, f"Negative template atom index {template_idx}"

    aligned = align_from_atom_pairs(source_mol, template_mol, atom_pairs)
    assert aligned is source_mol  # Should modify in place

    # Check that the alignment is correct
    _assert_coords_close(source_mol, template_mol, atom_pairs)

    # Test reverse alignment
    source_mol = Chem.MolFromSmiles(source_smiles)
    template_mol = Chem.MolFromSmiles(template_smiles)

    reversed_pairs = [(b, a) for a, b in atom_pairs]
    aligned = align_from_atom_pairs(template_mol, source_mol, reversed_pairs)
    _assert_coords_close(template_mol, source_mol, reversed_pairs)


def test_align_from_atom_pairs_invalid_indices():
    """Manual alignment with invalid indices."""
    template_mol = Chem.MolFromSmiles("CCC")  # 3 atoms
    source_mol = Chem.MolFromSmiles("CC")  # 2 atoms
    rdDepictor.Compute2DCoords(template_mol)
    rdDepictor.Compute2DCoords(source_mol)

    # Test template index out of bounds
    with pytest.raises(AssertionError, match="Template atom index .* out of range"):
        align_from_atom_pairs(source_mol, template_mol, [(0, template_mol.GetNumAtoms())])

    # Test molecule index out of bounds
    with pytest.raises(AssertionError, match="Source molecule atom index .* out of range"):
        align_from_atom_pairs(source_mol, template_mol, [(source_mol.GetNumAtoms(), 0)])

    # Test negative indices
    with pytest.raises(AssertionError, match="Negative .* atom index"):
        align_from_atom_pairs(source_mol, template_mol, [(-1, 0)])


@pytest.mark.parametrize(
    "source_smiles,template_smiles,index_mapping,atom_pairs",
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
def test_align_from_indices_valid(
    source_smiles, template_smiles, index_mapping, atom_pairs
):
    """Test align_from_indices"""
    source_mol = Chem.MolFromSmiles(source_smiles)
    template_mol = Chem.MolFromSmiles(template_smiles)

    assert source_mol, f"Failed to parse source molecule SMILES: {source_smiles}"
    assert template_mol, f"Failed to parse template SMILES: {template_smiles}"

    aligned = align_from_indices(source_mol, template_mol, index_mapping)
    assert aligned is source_mol  # Should modify in place

    # Check that the alignment is correct
    _assert_coords_close(source_mol, template_mol, atom_pairs)


@pytest.mark.parametrize(
    "source_smiles,template_smiles,index_mapping,error_type,error_match",
    [
        (
            "CC",
            "CC",
            [-1, 2],
            AssertionError,
            r"Template atom index 2 out of range \(max 1\)",
        ),  # Invalid template index
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
            "Length of index_mapping must equal",
        ),  # Wrong length list
    ],
    ids=["invalid_index", "invalid_type", "wrong_length"],
)
def test_align_from_indices_invalid(
    source_smiles, template_smiles, index_mapping, error_type, error_match
):
    """Test align_from_indices with invalid inputs"""
    source_mol = Chem.MolFromSmiles(source_smiles)
    template_mol = Chem.MolFromSmiles(template_smiles)

    assert source_mol, f"Failed to parse source molecule SMILES: {source_smiles}"
    assert template_mol, f"Failed to parse template SMILES: {template_smiles}"

    with pytest.raises(error_type, match=error_match):
        align_from_indices(source_mol, template_mol, index_mapping)


@pytest.mark.parametrize(
    "source_smiles,template_smiles,source_map_smiles,template_map_smiles,atom_pairs",
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
            "CC",  # 3-atom to 2-atom
            "[CH3:1][CH2:2][CH3:3]",
            "[CH3:1][CH3:2]",
            [(0, 0), (1, 1)],
        ),
        (
            "c1ccccc1",
            "c1ccccc1",  # Benzene ring
            "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1",
            "[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1",
            [(i, i) for i in range(6)],
        ),
    ],
    ids=["small_to_small", "large_to_small", "ring_to_ring"],
)
def test_align_from_mapids_valid(
    source_smiles, template_smiles, source_map_smiles, template_map_smiles, atom_pairs
):
    """Test align_from_mapids with valid inputs"""
    source_mol = Chem.MolFromSmiles(source_smiles)
    template_mol = Chem.MolFromSmiles(template_smiles)
    source_map = Chem.MolFromSmiles(source_map_smiles)
    template_map = Chem.MolFromSmiles(template_map_smiles)

    assert source_mol, f"Failed to parse source molecule SMILES: {source_smiles}"
    assert template_mol, f"Failed to parse template SMILES: {template_smiles}"
    assert source_map, f"Failed to parse source map SMILES: {source_map_smiles}"
    assert template_map, f"Failed to parse template map SMILES: {template_map_smiles}"

    aligned = align_from_mapids(source_mol, template_mol, source_map, template_map)
    assert aligned is source_mol  # Should modify in place

    # Check that the alignment is correct
    _assert_coords_close(source_mol, template_mol, atom_pairs)


@pytest.mark.parametrize(
    "source_smiles,template_smiles,source_map_smiles,template_map_smiles,error_type,error_match",
    [
        (
            "CC",
            "CC",
            "CC",  # No atom maps
            "[CH3:1][CH3:2]",
            ValueError,
            "No atom maps found in source molecule",
        ),
        (
            "CC",
            "CC",
            "[CH3:1][CH3:2]",
            "CC",  # No atom maps in template
            ValueError,
            "No atom maps found in template molecule",
        ),
        (
            "CC",
            "CC",
            "[CH3:1][CH3:3]",  # Non-matching map IDs
            "C[CH3:2]",  # Only has map ID 2, no overlap with source map IDs 1,3
            ValueError,
            "No matching atom map IDs found",
        ),
    ],
    ids=["no_source_maps", "no_template_maps", "non_matching_maps"],
)
def test_align_from_mapids_invalid(
    source_smiles,
    template_smiles,
    source_map_smiles,
    template_map_smiles,
    error_type,
    error_match,
):
    """Test align_from_mapids with invalid inputs"""
    source_mol = Chem.MolFromSmiles(source_smiles)
    template_mol = Chem.MolFromSmiles(template_smiles)
    source_map = Chem.MolFromSmiles(source_map_smiles)
    template_map = Chem.MolFromSmiles(template_map_smiles)

    assert source_mol, f"Failed to parse source molecule SMILES: {source_smiles}"
    assert template_mol, f"Failed to parse template SMILES: {template_smiles}"
    assert source_map, f"Failed to parse source map SMILES: {source_map_smiles}"
    assert template_map, f"Failed to parse template map SMILES: {template_map_smiles}"

    with pytest.raises(error_type, match=error_match):
        align_from_mapids(source_mol, template_mol, source_map, template_map)


def GetCoords(mol: Chem.Mol, i: int):
    c = mol.GetConformer(0).GetAtomPosition(i)
    return c.x, c.y, c.z


def _assert_coords_close(source_mol: Chem.Mol, template_mol: Chem.Mol, atom_pairs: List[Tuple[int, int]]):
    """Helper function to check if coordinates are aligned correctly.
    
    Args:
        source_mol: The source molecule to check
        template_mol: The template molecule to compare against
        atom_pairs: List of (source_idx, template_idx) pairs to compare
    """
    source_conf = source_mol.GetConformer()
    template_conf = template_mol.GetConformer()

    for source_idx, template_idx in atom_pairs:
        source_pos = source_conf.GetAtomPosition(source_idx)
        template_pos = template_conf.GetAtomPosition(template_idx)
        assert abs(source_pos.x - template_pos.x) < 0.1
        assert abs(source_pos.y - template_pos.y) < 0.1
        assert abs(source_pos.z - template_pos.z) < 0.1


def test_auto_align_molecules():
    """Test auto_align_molecules with a list of molecules."""
    molecules = [
        Chem.MolFromSmiles("CC"),  # ethane
        Chem.MolFromSmiles("CCC"),  # propane
        Chem.MolFromSmiles("CCCC"),  # butane
    ]
    assert all(molecules), "Failed to parse test molecules"

    # Test with default template (first molecule)
    aligned = auto_align_molecules(molecules)
    assert aligned == molecules  # Should modify in place and return input list

    # Test with explicit template via hint
    hint = Alignment.from_mcs(molecules[1], molecules[0])  # align propane to ethane
    aligned = auto_align_molecules(molecules, hints=[hint])
    assert aligned == molecules


def test_auto_align_molecules_with_hints():
    """Test auto_align_molecules with alignment hints."""

    # Create a few test molecules
    mol1 = Chem.MolFromSmiles("CCO")  # ethanol
    mol2 = Chem.MolFromSmiles("CCCO")  # propanol
    mol3 = Chem.MolFromSmiles("CCCCO")  # butanol

    mols = [mol1, mol2, mol3]

    # Generate 2D coords for first molecule
    rdDepictor.Compute2DCoords(mol1)

    # Create hint alignment between mol1 and mol2, different than the auto alignment
    atom_pairs = [(0, 0), (1, 1)]
    hint_alignment = Alignment(
        atom_pairs=atom_pairs,
        score=2.0,
        source_mol=mol1,
        template_mol=mol2
    )

    # Align all molecules
    aligned = auto_align_molecules(mols, hints=[hint_alignment])

    # Should return same number of molecules
    assert len(aligned) == len(mols)

    # Should modify in place
    assert aligned[0] is mol1
    assert aligned[1] is mol2
    assert aligned[2] is mol3

    # Check that mol1 and mol2 are aligned as specified by the hint
    # mol1 -> mol2
    _assert_coords_close(mol1, mol2, atom_pairs)
    # mol2 -> mol3
    _assert_coords_close(mol2, mol3, [(0, 1), (1, 2), (2, 3), (3, 4)])


def test_auto_align_molecules_no_align():
    """Test that molecules have valid coordinates even when not aligned."""
    # Create two simple molecules with no MCS alignment
    mol1 = Chem.MolFromSmiles("CC")
    mol2 = Chem.MolFromSmiles("N#N")
    
    # Ensure they start with no conformers
    assert mol1.GetNumConformers() == 0
    assert mol2.GetNumConformers() == 0
    
    mols = auto_align_molecules([mol1, mol2])
    
    # Verify both molecules have conformers
    assert mols[0].GetNumConformers() > 0
    assert mols[1].GetNumConformers() > 0
    