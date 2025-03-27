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

# Test molecules
TEMPLATE_SMILES = 'c1ccccc1'  # benzene
MOL_SMILES = 'Cc1ccccc1'  # toluene
COMPLEX_TEMPLATE = 'c1nccc2n1ccc2'  # pyrrolopyridine
COMPLEX_MOL = 'OCCc1ccn2cnccc12'  # substituted pyrrolopyridine

def test_align_to_template():
    """Test basic automatic template alignment."""
    template = Chem.MolFromSmiles(TEMPLATE_SMILES)
    mol = Chem.MolFromSmiles(MOL_SMILES)
    
    # Test basic alignment
    aligned = align_to_template(mol, template)
    assert aligned.GetNumConformers() == 1
    assert aligned is mol  # Should modify in place
    
    # Test with more complex molecules
    template = Chem.MolFromSmiles(COMPLEX_TEMPLATE)
    mol = Chem.MolFromSmiles(COMPLEX_MOL)
    aligned = align_to_template(mol, template)
    assert aligned.GetNumConformers() == 1
    
    # Test with no matching substructure
    template = Chem.MolFromSmiles('CCO')  # ethanol
    mol = Chem.MolFromSmiles('c1ccccc1')  # benzene
    with pytest.raises(ValueError, match="Substructure match with reference not found"):
        align_to_template(mol, template)

def test_align_to_template_manual():
    """Test alignment using explicit atom pairs."""
    # Create fresh molecules each time
    template = Chem.MolFromSmiles('c1ccccc1')  # benzene
    mol = Chem.MolFromSmiles('Cc1ccccc1')  # toluene
    
    # Validate molecule creation
    assert template is not None, "Failed to create template molecule"
    assert mol is not None, "Failed to create test molecule"
    
    # Generate initial 2D coords
    rdDepictor.Compute2DCoords(template)
    rdDepictor.Compute2DCoords(mol)

    # Test with minimal valid atom pairs
    pairs = [(1,0), (2,1), (3,2)]
    aligned = align_to_template_manual(mol, template, pairs)
    assert aligned.GetNumConformers() == 1
    
    # Test with even fewer pairs
    pairs = [(1,0), (2,1)]
    aligned = align_to_template_manual(mol, template, pairs)
    assert aligned.GetNumConformers() == 1
    
    # Test with single pair
    pairs = [(1,0)]
    aligned = align_to_template_manual(mol, template, pairs)
    assert aligned.GetNumConformers() == 1
    
    # Test with empty pairs list
    aligned = align_to_template_manual(mol, template, [])
    assert aligned.GetNumConformers() == 1
    
    # Test with invalid indices should raise
    with pytest.raises(RuntimeError):  # RDKit raises RuntimeError for invalid indices
        align_to_template_manual(mol, template, [(mol.GetNumAtoms(),0)])

def test_align_to_template_with_indices():
    """Test alignment using atom index mappings."""
    template = Chem.MolFromSmiles(TEMPLATE_SMILES)  # 6 atoms
    mol = Chem.MolFromSmiles(MOL_SMILES)  # 7 atoms (including methyl)
    
    # Test with list mapping
    # Map 5 ring carbons of toluene to benzene (methyl and one carbon unmapped)
    mapping = [-1] * 6  # Initialize all unmapped
    for i in range(5):  # Map 5 atoms to ensure valid indices
        mapping[i] = i + 1  # Map template atom i to mol atom i+1
    aligned = align_to_template_with_indices(mol, template, mapping)
    assert aligned.GetNumConformers() == 1
    
    # Test with dict mapping
    mapping = {i+1: i for i in range(5)}  # Same mapping as above
    aligned = align_to_template_with_indices(mol, template, mapping)
    assert aligned.GetNumConformers() == 1
    
    # Test with empty mapping
    mapping = {i: -1 for i in range(6)}  # All unmapped
    aligned = align_to_template_with_indices(mol, template, mapping)
    assert aligned.GetNumConformers() == 1
    
    # Test with wrong list length
    with pytest.raises(AssertionError):
        align_to_template_with_indices(mol, template, [-1])
    
    # Test with invalid indices
    with pytest.raises(RuntimeError):
        align_to_template_with_indices(mol, template, [-2] * template.GetNumAtoms())

def test_align_to_template_by_mapids():
    """Test alignment using atom map IDs."""
    # Create molecules with atom mapping
    # Map 5 ring carbons, leaving one unmapped for valid indices
    template_map = Chem.MolFromSmiles('c1c([*:2])c([*:3])c([*:4])c([*:5])c1')
    mol_map = Chem.MolFromSmiles('C([*:1])c1c([*:2])c([*:3])c([*:4])c([*:5])c1')
    template = Chem.MolFromSmiles(TEMPLATE_SMILES)
    mol = Chem.MolFromSmiles(MOL_SMILES)
    
    # Test alignment with map IDs
    aligned = align_to_template_by_mapids(mol, template, mol_map, template_map)
    assert aligned.GetNumConformers() == 1
    
    # Test with no matching map IDs
    template_map = Chem.MolFromSmiles('c1c([*:7])c([*:8])c([*:9])c([*:10])c1')
    aligned = align_to_template_by_mapids(mol, template, mol_map, template_map)
    assert aligned.GetNumConformers() == 1  # Should still generate coords
    
    # Test with invalid molecules should raise
    with pytest.raises(ValueError):
        align_to_template_by_mapids(mol, template, 
                                  Chem.MolFromSmiles('invalid'), 
                                  template_map) 