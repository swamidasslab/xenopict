"""Declarative API for xenopict."""

from typing import List, Union
from ..types.base import XenopictSpec
from .. import Xenopict

def from_spec(spec: Union[dict, XenopictSpec]) -> List[str]:
    """Convert a specification into a list of SVG images.
    
    Args:
        spec: Either a dict matching the XenopictSpec schema or a XenopictSpec instance
        
    Returns:
        List of SVG strings, one for each molecule in the specification
    """
    if isinstance(spec, dict):
        spec = XenopictSpec(**spec)
    
    # Convert molecule field to list if it's a single molecule
    molecules = spec.molecule if isinstance(spec.molecule, list) else [spec.molecule]
    
    # Create SVG for each molecule
    images = []
    for mol in molecules:
        # Create Xenopict object from SMILES
        xeno = Xenopict(mol.smiles)
        # Convert to SVG and append to list
        images.append(xeno.to_svg())
    
    return images

def from_json(json_str: str) -> List[str]:
    """Convert a JSON string specification into a list of SVG images.
    
    Args:
        json_str: JSON string matching the XenopictSpec schema
        
    Returns:
        List of SVG strings, one for each molecule in the specification
    """
    import json
    spec = json.loads(json_str)
    return from_spec(spec)

def from_file(path: str) -> List[str]:
    """Load a specification from a JSON file and convert it into a list of SVG images.
    
    Args:
        path: Path to JSON file containing a specification matching the XenopictSpec schema
        
    Returns:
        List of SVG strings, one for each molecule in the specification
    """
    with open(path) as f:
        return from_json(f.read()) 