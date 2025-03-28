"""Declarative API for xenopict.

This module provides a high-level interface for creating molecule visualizations
using JSON-like specifications. It supports both single and multiple molecule
visualizations with automatic alignment.
"""

from typing import List, Union, Dict, Any, TypeVar, cast, overload
from pathlib import Path
from .types import XenopictSpec, MoleculeSpec
from .. import Xenopict
from ..alignment import auto_align_molecules, _ensure_coords
from rdkit.Chem import MolFromSmiles, Mol  # type: ignore
import json

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
        
        >>> # From a JSON string
        >>> json_str = '''
        ... {
        ...     "molecules": [
        ...         {
        ...             "smiles": "CCO"
        ...         },
        ...         {
        ...             "smiles": "CCCO"
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
    molecules = xenopict_spec.molecules if isinstance(xenopict_spec.molecules, list) else [xenopict_spec.molecules]
    
    # Create RDKit molecules
    rdkit_mols = [_process_smiles(mol.smiles) for mol in molecules]
    
    # Generate 2D coordinates for all molecules
    for mol in rdkit_mols:
        _ensure_coords(mol)
    
    # Align molecules if more than one and alignment not explicitly disabled
    if len(rdkit_mols) > 1 and xenopict_spec.align:
        auto_align_molecules(rdkit_mols)
    
    # Create Xenopict objects after alignment
    return [Xenopict(mol) for mol in rdkit_mols] 