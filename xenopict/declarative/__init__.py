"""Declarative API for xenopict."""

from typing import List, Union
from ..types.base import XenopictSpec
from .. import Xenopict

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
    
    # Create Xenopict object for each molecule
    return [Xenopict(mol.smiles) for mol in molecules]

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