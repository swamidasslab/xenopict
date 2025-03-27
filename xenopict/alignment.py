"""Functions for aligning 2D molecular depictions using RDKit.

This module provides several methods for aligning molecules to template structures:
- Basic alignment using automatic substructure matching
- Manual alignment using explicit atom pairs
- Alignment using atom indices
- Alignment using atom map IDs
"""

from rdkit import Chem  # type: ignore
from rdkit.Chem import rdDepictor  # type: ignore
from typing import List, Dict, Tuple, Union


def _ensure_coords(mol: Chem.Mol) -> Chem.Mol:
    """Ensure that the molecule has 2D coordinates.
    
    Args:
        mol: The molecule to check/compute coordinates for
        
    Returns:
        The input molecule (modified in place with 2D coords if needed)
    """
    if mol.GetNumConformers() == 0:
        rdDepictor.Compute2DCoords(mol)
    return mol


def align_to_template(mol: Chem.Mol, template: Chem.Mol) -> Chem.Mol:
    """Align a molecule to a template molecule using automatic substructure matching.

    This function attempts to align the input molecule to match the 2D depiction
    of the template molecule. The alignment is done automatically using RDKit's
    substructure matching.

    Args:
        mol: The molecule to align (modified in place)
        template: The template molecule to align to
        
    Returns:
        The aligned molecule (same object as input mol)
    """
    _ensure_coords(template)
    rdDepictor.GenerateDepictionMatching2DStructure(mol, template)
    return mol


def align_to_template_manual(
    mol: Chem.Mol, 
    template: Chem.Mol, 
    aligned_atoms: List[Tuple[int, int]]
) -> Chem.Mol:
    """Align a molecule to a template molecule using explicit atom pairs.

    This function aligns the input molecule to match the 2D depiction of the template
    molecule using explicitly specified atom pairs.

    Args:
        mol: The molecule to align (modified in place)
        template: The template molecule to align to
        aligned_atoms: List of (mol_atom_idx, template_atom_idx) pairs specifying 
                      which atoms should be aligned to each other
        
    Returns:
        The aligned molecule (same object as input mol)
    """
    _ensure_coords(template)
    if len(aligned_atoms):  # empty list causes segfault
        rdDepictor.GenerateDepictionMatching2DStructure(mol, template, aligned_atoms)
    else:
        _ensure_coords(mol)
    return mol


def align_to_template_with_indices(
    mol: Chem.Mol, 
    template: Chem.Mol, 
    mol_to_template: Union[List[int], Dict[int, int]]
) -> Chem.Mol:
    """Align a molecule to a template molecule using atom index mappings.

    This function aligns the input molecule to match the 2D depiction of the template
    molecule using a mapping between atom indices.

    Args:
        mol: The molecule to align (modified in place)
        template: The template molecule to align to
        mol_to_template: Either:
            - A list where index i contains the template atom index that mol atom i maps to
              (-1 for unmapped atoms). Length must match number of atoms in template.
            - A dict mapping mol atom indices to template atom indices
              (omit or use -1 for unmapped atoms)
        
    Returns:
        The aligned molecule (same object as input mol)
        
    Raises:
        ValueError: If mol_to_template is not a list or dict
        AssertionError: If mol_to_template is a list with wrong length
    """
    _ensure_coords(template)

    if isinstance(mol_to_template, list):
        assert len(mol_to_template) == template.GetNumAtoms(), (
            f"Length of mol_to_template must equal number of template atoms: "
            f"{len(mol_to_template)} != {template.GetNumAtoms()}"
        )

        # Convert list format to explicit atom pairs, skipping unmapped atoms (-1)
        aligned_atoms = [
            (mol_idx, template_idx)
            for mol_idx, template_idx in enumerate(mol_to_template)
            if template_idx >= 0
        ]

    elif isinstance(mol_to_template, dict):
        # Convert dict format to explicit atom pairs, skipping unmapped atoms (-1)
        aligned_atoms = [
            (mol_idx, template_idx)
            for mol_idx, template_idx in mol_to_template.items()
            if template_idx >= 0
        ]
    else:
        raise ValueError(
            f"mol_to_template must be a list or dictionary, got {type(mol_to_template)}"
        )

    if len(aligned_atoms):  # empty list causes segfault
        rdDepictor.GenerateDepictionMatching2DStructure(mol, template, aligned_atoms)
    else:
        _ensure_coords(mol)
    return mol


def _extract_map_ids(mol: Chem.Mol) -> Dict[int, int]:
    """Extract atom map IDs from a molecule.
    
    Args:
        mol: The molecule to extract map IDs from
        
    Returns:
        Dictionary mapping atom map IDs to atom indices
    """
    return {
        atom.GetAtomMapNum(): atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomMapNum() != 0
    }


def _get_matched_mapids(mol: Chem.Mol, mol_map: Chem.Mol) -> Dict[int, int]:
    """Get mapping between atom map IDs and atom indices for matched atoms.

    Args:
        mol: The molecule to match against
        mol_map: The molecule containing the atom mapping IDs to match

    Returns:
        Dictionary mapping atom map IDs to matched atom indices in mol
    """
    # Get map ID to index mapping from the mapping molecule
    mapid2idx = _extract_map_ids(mol_map)
    
    # Find substructure match between molecules
    match = mol.GetSubstructMatch(mol_map)
    
    # Map the map IDs to the matched atom indices
    return {
        map_id: match[idx] 
        for map_id, idx in mapid2idx.items() 
        if match[idx] >= 0
    }

  
def align_to_template_by_mapids(
    mol: Chem.Mol, 
    template: Chem.Mol, 
    mol_map: Chem.Mol, 
    template_map: Chem.Mol
) -> Chem.Mol:
    """Align a molecule to a template using atom map IDs.

    This function aligns the input molecule to match the 2D depiction of the template
    molecule using atom map IDs to identify corresponding atoms. The mapping is done
    through two auxiliary molecules (mol_map and template_map) that contain the map IDs.

    Args:
        mol: The molecule to align (modified in place)
        template: The template molecule to align to
        mol_map: Molecule with map IDs matching atoms in mol
        template_map: Molecule with map IDs matching atoms in template

    Returns:
        The aligned molecule (same object as input mol)
    """
    _ensure_coords(template)

    # Get map ID to atom index mappings for both molecules
    mapid2idx_mol = _get_matched_mapids(mol, mol_map)
    mapid2idx_template = _get_matched_mapids(template, template_map)

    # Find map IDs present in both molecules
    overlapping_mapids = set(mapid2idx_mol) & set(mapid2idx_template)

    # Create list of atom pairs for alignment
    aligned_atoms = [
        (mapid2idx_mol[mapid], mapid2idx_template[mapid])
        for mapid in overlapping_mapids
    ]

    if len(aligned_atoms):  # empty list causes segfault
        rdDepictor.GenerateDepictionMatching2DStructure(mol, template, aligned_atoms)
    else:
        _ensure_coords(mol)
    return mol
