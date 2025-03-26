"""Declarative API for xenopict."""

from typing import List, Union, Sequence, cast, Tuple
from ..types.base import XenopictSpec, MoleculeSpec, ShadeSpec, CircleSpec
from .. import Xenopict
from rdkit.Chem import MolFromSmarts, MolFromSmiles  # type: ignore

def _apply_shading(xeno: Xenopict, shade_spec: ShadeSpec) -> None:
    """Apply shading to a Xenopict object based on a ShadeSpec."""
    if shade_spec.type == "atom":
        # For atom shading, targets should be a list of atom indices
        atom_indices = cast(List[int], shade_spec.targets)
        values = [float(shade_spec.value)] * len(atom_indices)
        xeno.shade(atom_shading=values)  # type: ignore
    elif shade_spec.type == "bond":
        # For bond shading, targets should be list of [atom1, atom2] pairs
        bond_pairs = cast(List[List[int]], shade_spec.targets)
        atom1 = [pair[0] for pair in bond_pairs]
        atom2 = [pair[1] for pair in bond_pairs]
        values = [float(shade_spec.value)] * len(bond_pairs)
        xeno.shade(bond_shading=(atom1, atom2, values))
    elif shade_spec.type == "substructure":
        # For substructure, targets should be a SMARTS pattern
        smarts = cast(str, shade_spec.targets)
        patt = MolFromSmarts(smarts)
        if patt is not None:
            hit_ats = xeno.mol.GetSubstructMatch(patt)
            if hit_ats:
                xeno.shade_substructure([list(hit_ats)], [float(shade_spec.value)])

def _apply_circles(xeno: Xenopict, circle_spec: CircleSpec) -> None:
    """Apply circle annotations to a Xenopict object based on a CircleSpec."""
    if circle_spec.type == "atom":
        # For atoms, targets should be a list of atom indices
        atom_indices = cast(List[int], circle_spec.targets)
        xeno.mark_substructure(atom_indices)
    elif circle_spec.type == "bond":
        # For bonds, targets should be list of [atom1, atom2] pairs
        bond_pairs = cast(List[List[int]], circle_spec.targets)
        xeno.mark_substructure([], bond_pairs)

def _process_molecule(mol_spec: MoleculeSpec) -> Xenopict:
    """Process a single molecule specification into a Xenopict object."""
    # Validate that either SMILES or SMARTS is provided
    if mol_spec.smiles is None and mol_spec.smarts is None:
        raise ValueError("Either 'smiles' or 'smarts' must be provided")
    
    # Convert to RDKit mol using appropriate function
    if mol_spec.smiles is not None:
        mol = MolFromSmiles(mol_spec.smiles)
        if mol is None:
            raise ValueError(f"Invalid SMILES: {mol_spec.smiles}")
    else:
        # Must be SMARTS since we validated above
        mol = MolFromSmarts(mol_spec.smarts)  # type: ignore
        if mol is None:
            raise ValueError(f"Invalid SMARTS pattern: {mol_spec.smarts}")  # type: ignore
    
    # Create Xenopict object from RDKit mol
    xeno = Xenopict(mol)
    
    # Apply backbone color if specified
    if mol_spec.backbone_color:
        xeno.set_backbone_color(mol_spec.backbone_color)
    
    # Apply shading if specified
    if mol_spec.shading:
        for shade in mol_spec.shading:
            _apply_shading(xeno, shade)
    
    # Apply circles if specified
    if mol_spec.circles:
        for circle in mol_spec.circles:
            _apply_circles(xeno, circle)
    
    return xeno

def from_spec(spec: Union[dict, XenopictSpec]) -> List[Xenopict]:
    """Convert a specification into a list of Xenopict objects.
    
    Args:
        spec: Either a dict matching the XenopictSpec schema or a XenopictSpec instance
        
    Returns:
        List of Xenopict objects, one for each molecule in the specification
    """
    if isinstance(spec, dict):
        spec = XenopictSpec(**spec)
    
    # Convert molecule field to list if it's a single molecule
    molecules = spec.molecule if isinstance(spec.molecule, list) else [spec.molecule]
    
    # Process each molecule
    return [_process_molecule(mol) for mol in molecules]

def from_json(json_str: str) -> List[Xenopict]:
    """Convert a JSON string specification into a list of Xenopict objects.
    
    Args:
        json_str: JSON string matching the XenopictSpec schema
        
    Returns:
        List of Xenopict objects, one for each molecule in the specification
    """
    import json
    spec = json.loads(json_str)
    return from_spec(spec)

def from_file(path: str) -> List[Xenopict]:
    """Load a specification from a JSON file and convert it into a list of Xenopict objects.
    
    Args:
        path: Path to JSON file containing a specification matching the XenopictSpec schema
        
    Returns:
        List of Xenopict objects, one for each molecule in the specification
    """
    with open(path) as f:
        return from_json(f.read()) 