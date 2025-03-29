"""Declarative API for xenopict.

This module provides a high-level interface for creating molecule visualizations
using JSON-like specifications. It supports both single and multiple molecule
visualizations with automatic alignment.
"""

# NOTE: This module requires Xenopict to be imported first in xenopict/__init__.py
# ruff: noqa: I001

import json
from pathlib import Path
from typing import Any, Dict, List, TypeVar, Union, overload

from rdkit.Chem import Mol, MolFromSmiles  # type: ignore

from .. import Xenopict
from ..alignment import _ensure_coords, auto_align_molecules
from .types import MarkSpec, MoleculeSpec, XenopictSpec

T = TypeVar("T", bound=Union[Dict[str, Any], XenopictSpec])


def _process_smiles(smiles: str) -> Mol:
    """Convert SMILES to RDKit molecule.

    Args:
        smiles: SMILES string to convert

    Returns:
        RDKit molecule object

    Raises:
        ValueError: If the SMILES string is invalid
    """
    mol = MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return mol


def _apply_marks(xenopict: Xenopict, mark_spec: MarkSpec) -> None:
    """Apply marking specifications to a Xenopict object.

    Args:
        xenopict: The Xenopict object to mark
        mark_spec: The marking specifications to apply
    """
    if mark_spec.substructure_atoms is not None:
        xenopict.mark_substructure(mark_spec.substructure_atoms, mark_spec.substructure_bonds)
    if mark_spec.atoms is not None:
        xenopict.mark_atoms(mark_spec.atoms)


def _apply_style(xenopict: Xenopict, mol_spec: MoleculeSpec, xenopict_spec: XenopictSpec) -> None:
    """Apply style specifications to a Xenopict object.

    Args:
        xenopict: The Xenopict object to style
        mol_spec: The molecule specification containing style options
        xenopict_spec: The parent specification containing default style options
    """
    # Apply color - molecule specific overrides default
    color = mol_spec.color if mol_spec.color is not None else xenopict_spec.color
    if color is not None:
        xenopict.set_backbone_color(color)

    # Apply halo - molecule specific overrides default
    # If mol_spec.halo is None, use xenopict_spec.halo (which defaults to True)
    use_halo = mol_spec.halo if mol_spec.halo is not None else xenopict_spec.halo
    if use_halo:
        xenopict.halo()  # Add halo elements


@overload
def parse(spec: Union[Dict[str, Any], XenopictSpec]) -> List[Xenopict]: ...


@overload
def parse(spec: str) -> List[Xenopict]: ...


@overload
def parse(spec: Path) -> List[Xenopict]: ...


def parse(spec: Union[Dict[str, Any], XenopictSpec, str, Path]) -> List[Xenopict]:
    """Parse a specification into a list of Xenopict objects.

    This function serves as the main entry point for the declarative API. It accepts
    specifications in multiple formats and returns a list of Xenopict objects.

    Args:
        spec: The specification for creating molecules. Can be one of:
            - A dict matching the XenopictSpec schema
            - A XenopictSpec instance
            - A JSON string matching the XenopictSpec schema
            - A Path to a JSON file containing a XenopictSpec

    Returns:
        List of Xenopict objects, one for each molecule in the specification

    Raises:
        ValueError: If the input is invalid or required fields are missing

    Examples:
        >>> # From a dict
        >>> spec = {
        ...     "molecules": {
        ...         "smiles": "CCO"
        ...     }
        ... }
        >>> xenopicts = parse(spec)
        >>> len(xenopicts)
        1

        >>> # From a JSON string with marking
        >>> json_str = '''
        ... {
        ...     "molecules": [
        ...         {
        ...             "smiles": "CCO",
        ...             "mark": {
        ...                 "atoms": [0, 1]
        ...             }
        ...         },
        ...         {
        ...             "smiles": "CCCO",
        ...             "mark": {
        ...                 "substructure_atoms": [0, 1, 2]
        ...             }
        ...         }
        ...     ]
        ... }
        ... '''
        >>> xenopicts = parse(json_str)
        >>> len(xenopicts)
        2

        >>> # From a file path
        >>> from pathlib import Path
        >>> path = Path("molecules.json")  # doctest: +SKIP
        >>> xenopicts = parse(path)  # doctest: +SKIP
    """
    # Convert input to XenopictSpec
    if isinstance(spec, (dict, XenopictSpec)):
        xenopict_spec = spec if isinstance(spec, XenopictSpec) else XenopictSpec(**spec)
    elif isinstance(spec, str):
        try:
            xenopict_spec = XenopictSpec(**json.loads(spec))
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON string: {str(e)}")
    elif isinstance(spec, Path):
        try:
            with open(spec) as f:
                xenopict_spec = XenopictSpec(**json.loads(f.read()))
        except FileNotFoundError:
            raise ValueError(f"File not found: {spec}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Invalid JSON in file {spec}: {str(e)}")
        except IOError as e:
            raise ValueError(f"Error reading file {spec}: {str(e)}")
    else:
        raise TypeError(f"Unsupported input type: {type(spec)}")

    # Convert single molecule to list for consistent handling
    molecules = (
        xenopict_spec.molecules
        if isinstance(xenopict_spec.molecules, list)
        else [xenopict_spec.molecules]
    )

    # Create RDKit molecules
    rdkit_mols = [_process_smiles(mol.smiles) for mol in molecules]

    # Generate 2D coordinates for all molecules
    for mol in rdkit_mols:
        _ensure_coords(mol)

    # Align molecules if more than one and alignment not explicitly disabled
    if len(rdkit_mols) > 1 and xenopict_spec.align:
        auto_align_molecules(rdkit_mols)

    # Create Xenopict objects after alignment
    xenopicts = [Xenopict(mol) for mol in rdkit_mols]

    # Apply marking specifications and styles if present
    for xenopict, mol_spec in zip(xenopicts, molecules):
        if mol_spec.mark is not None:
            _apply_marks(xenopict, mol_spec.mark)
        _apply_style(
            xenopict, mol_spec, xenopict_spec
        )  # Always call _apply_style to handle defaults

    return xenopicts
