"""Tests for molecular alignment functionality."""

from rdkit import Chem
import pytest
from rdkit.Chem import rdDepictor

from xenopict.alignment import (
    align_to_template,
    align_to_template_manual,
    align_to_template_with_indices,
    align_to_template_by_mapids,
)

# Test molecules and their descriptions for parameterized tests
ALIGNMENT_TEST_CASES = [
    ('c1ccccc1', 'Cc1ccccc1', 'benzene->toluene'),
    ('c1nccc2n1ccc2', 'OCCc1ccn2cnccc12', 'pyrrolopyridine->subst'),
]

MANUAL_PAIRS_TEST_CASES = [
    ([(1,0), (2,1), (3,2)], "3 pairs"),
    ([(1,0), (2,1)], "2 pairs"),
    ([(1,0)], "1 pair"),
    ([], "no pairs")
]

INDEX_MAPPING_TEST_CASES = [
    # (mapping, description)
    ([-1, 1, 2, 3, 4, 5], "list-5atoms"),
    ({1: 0, 2: 1, 3: 2, 4: 3, 5: 4}, "dict-5atoms"),
    ({i: -1 for i in range(6)}, "empty-dict")
]

@pytest.fixture
def basic_mols():
    """Basic benzene/toluene test pair with 2D coords."""
    template = Chem.MolFromSmiles('c1ccccc1')
    mol = Chem.MolFromSmiles('Cc1ccccc1')
    rdDepictor.Compute2DCoords(template)
    rdDepictor.Compute2DCoords(mol)
    return template, mol

@pytest.mark.parametrize("template_smiles,mol_smiles,description", ALIGNMENT_TEST_CASES)
def test_align_to_template(template_smiles, mol_smiles, description):
    """Basic template alignment."""
    template = Chem.MolFromSmiles(template_smiles)
    mol = Chem.MolFromSmiles(mol_smiles)
    
    aligned = align_to_template(mol, template)
    assert aligned.GetNumConformers() == 1
    assert aligned is mol  # Should modify in place

def test_align_to_template_no_match():
    """Fails when no substructure match."""
    template = Chem.MolFromSmiles('CCO')  # ethanol
    mol = Chem.MolFromSmiles('c1ccccc1')  # benzene
    with pytest.raises(ValueError, match="Substructure match with reference not found"):
        align_to_template(mol, template)

@pytest.mark.parametrize("pairs,description", MANUAL_PAIRS_TEST_CASES)
def test_align_to_template_manual_pairs(basic_mols, pairs, description):
    """Manual alignment with atom pairs."""
    template, mol = basic_mols
    aligned = align_to_template_manual(mol, template, pairs)
    assert aligned.GetNumConformers() == 1

def test_align_to_template_manual_invalid_indices(basic_mols):
    """Manual alignment with invalid indices."""
    template, mol = basic_mols
    with pytest.raises(ValueError, match="Reference atom index in refMatchVect out of range"):
        align_to_template_manual(mol, template, [(mol.GetNumAtoms(), 0)])

@pytest.mark.parametrize("mapping,description", INDEX_MAPPING_TEST_CASES)
def test_align_to_template_with_indices_valid(basic_mols, mapping, description):
    """Index-based alignment."""
    template, mol = basic_mols
    aligned = align_to_template_with_indices(mol, template, mapping)
    assert aligned.GetNumConformers() == 1

def test_align_to_template_with_indices_invalid(basic_mols):
    """Index-based alignment with invalid inputs."""
    template, mol = basic_mols
    with pytest.raises(AssertionError):
        align_to_template_with_indices(mol, template, [-1])  # Wrong length
    with pytest.raises(ValueError, match="Reference atom index in refMatchVect out of range"):
        align_to_template_with_indices(mol, template, [template.GetNumAtoms()] * template.GetNumAtoms())  # Too large indices

@pytest.fixture
def mapped_mols():
    """Molecules with atom mapping."""
    template_map = Chem.MolFromSmiles('c1c([*:2])c([*:3])c([*:4])c([*:5])c1')
    mol_map = Chem.MolFromSmiles('C([*:1])c1c([*:2])c([*:3])c([*:4])c([*:5])c1')
    template = Chem.MolFromSmiles('c1ccccc1')
    mol = Chem.MolFromSmiles('Cc1ccccc1')
    rdDepictor.Compute2DCoords(template)
    rdDepictor.Compute2DCoords(mol)
    return template, mol, template_map, mol_map

def test_align_to_template_by_mapids_matching(mapped_mols):
    """Alignment with matching map IDs."""
    template, mol, template_map, mol_map = mapped_mols
    aligned = align_to_template_by_mapids(mol, template, mol_map, template_map)
    assert aligned.GetNumConformers() == 1

def test_align_to_template_by_mapids_no_match(mapped_mols):
    """Alignment with non-matching map IDs."""
    template, mol, _, mol_map = mapped_mols
    template_map = Chem.MolFromSmiles('c1c([*:7])c([*:8])c([*:9])c([*:10])c1')
    aligned = align_to_template_by_mapids(mol, template, mol_map, template_map)
    assert aligned.GetNumConformers() == 1

def test_align_to_template_by_mapids_invalid(mapped_mols):
    """Alignment with invalid mapped molecules."""
    template, mol, template_map, _ = mapped_mols
    with pytest.raises(ValueError, match="Invalid SMILES"):
        align_to_template_by_mapids(mol, template, 
                                  Chem.MolFromSmiles('invalid'), 
                                  template_map) 