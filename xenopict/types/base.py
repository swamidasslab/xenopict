"""Base type definitions for xenopict's declarative API."""

from typing import Dict, List, Optional, Union, Literal
from pydantic import BaseModel, Field, field_validator, confloat
from enum import Enum
import json

class AlignmentMethod(str, Enum):
    """Methods for aligning molecules."""
    AUTO = "auto"  # Let RDKit automatically align
    SUBSTRUCTURE = "substructure"  # Align based on common substructure
    MAPIDS = "mapids"  # Align based on atom map IDs in SMILES
    MANUAL = "manual"  # Manual alignment using atom indices

class AlignmentSpec(BaseModel):
    """Specification for molecule alignment."""
    method: AlignmentMethod
    reference_mol: str  # ID of reference molecule
    target_mol: str  # ID of molecule to align
    params: Optional[Dict] = Field(
        None,
        description="Method-specific parameters",
        examples=[
            # For SUBSTRUCTURE
            {"smarts": "c1ccccc1"},  
            # For MAPIDS - no params needed, uses map IDs in SMILES
            None,
            # For MANUAL
            {"atom_pairs": [[1, 1], [2, 2], [3, 3]]}  # List of [ref_atom_idx, target_atom_idx]
        ]
    )

class ShadeType(str, Enum):
    """Types of shading targets."""
    ATOM = "atom"
    BOND = "bond"
    SUBSTRUCTURE = "substructure"

class ShadeSpec(BaseModel):
    """Specification for shading parts of a molecule with a heatmap value."""
    type: ShadeType
    value: float = Field(
        ...,
        description="Value for heatmap shading (0 to 1 or -1 to 1)",
        ge=-1.0,
        le=1.0
    )
    targets: Union[List[int], List[List[int]], str] = Field(
        ...,
        description="Atoms indices, bond indices pairs, or SMARTS pattern",
        examples=[
            [1, 2, 3],  # Atom indices
            [[1, 2], [3, 4]],  # Bond indices
            "c1ccccc1"  # SMARTS pattern
        ]
    )
    zero_centered: bool = Field(
        False,
        description="If True, value range is -1 to 1, otherwise 0 to 1"
    )

class CircleSpec(BaseModel):
    """Specification for circling atoms or bonds."""
    type: Literal["atom", "bond"]
    targets: Union[List[int], List[List[int]]] = Field(
        ...,
        description="Atoms indices or bond indices pairs to circle",
        examples=[
            [1, 2, 3],  # Atom indices
            [[1, 2], [3, 4]]  # Bond indices
        ]
    )
    color: str = Field("#FF0000", pattern="^#[0-9a-fA-F]{6}$")
    width: float = Field(1.0, description="Line width of the circle")

class MoleculeSpec(BaseModel):
    """Specification for a single molecule."""
    id: Optional[str] = Field(None, description="Optional identifier for the molecule, required only when referenced by alignments or other features")
    smiles: Optional[str] = Field(None, description="SMILES string of the molecule")
    smarts: Optional[str] = Field(None, description="SMARTS pattern of the molecule")
    shading: Optional[List[ShadeSpec]] = None
    circles: Optional[List[CircleSpec]] = None
    backbone_color: Optional[str] = Field(
        None, 
        pattern="^#[0-9a-fA-F]{6}$",
        description="Color for the molecule's backbone"
    )

    @field_validator("smarts")
    @classmethod
    def validate_smarts_smiles_mutually_exclusive(cls, v: Optional[str], info) -> Optional[str]:
        """Validate that only one of SMILES or SMARTS is provided."""
        if v is not None and info.data.get("smiles") is not None:
            raise ValueError("Only one of 'smiles' or 'smarts' should be provided")
        return v

    @field_validator("smiles", "smarts")
    @classmethod
    def validate_smiles_or_smarts_required(cls, v: Optional[str], info) -> Optional[str]:
        """Validate that either SMILES or SMARTS is provided."""
        if info.field_name == "smarts" and v is None and info.data.get("smiles") is None:
            raise ValueError("Either 'smiles' or 'smarts' must be provided")
        return v

class XenopictSpec(BaseModel):
    """Root specification for xenopict visualizations."""
    molecule: Union[MoleculeSpec, List[MoleculeSpec]] = Field(
        ...,
        description="Single molecule or list of molecules to render"
    )
    alignments: Optional[List[AlignmentSpec]] = None

    @field_validator("alignments")
    @classmethod
    def validate_alignments(cls, v: Optional[List[AlignmentSpec]], info) -> Optional[List[AlignmentSpec]]:
        """Validate that referenced molecules have IDs when alignments are used."""
        if v is not None:
            molecules = info.data.get("molecule")
            if not isinstance(molecules, list):
                molecules = [molecules]
            
            # Create a set of available IDs
            molecule_ids = {m.id for m in molecules if m.id is not None}
            
            # Check each alignment
            for alignment in v:
                if alignment.reference_mol not in molecule_ids:
                    raise ValueError(f"Reference molecule '{alignment.reference_mol}' not found or missing ID")
                if alignment.target_mol not in molecule_ids:
                    raise ValueError(f"Target molecule '{alignment.target_mol}' not found or missing ID")
        return v

    class Config:
        """Pydantic model configuration with examples."""
        json_schema_extra = {
            "examples": [
                # Example 1: Single molecule
                {
                    "molecule": {
                        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
                    }
                },
                
                # Example 2: Multiple molecules with alignment
                {
                    "molecule": [
                        {
                            "id": "mol1",
                            "smiles": "[CH3:1][CH2:2][OH:3]"
                        },
                        {
                            "id": "mol2",
                            "smiles": "[CH3:1][CH2:2][NH2:3]"
                        }
                    ],
                    "alignments": [
                        {
                            "method": "mapids",
                            "reference_mol": "mol1",
                            "target_mol": "mol2"
                        }
                    ]
                },
                
                # Example 3: Single annotated molecule
                {
                    "molecule": {
                        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                        "shading": [
                            {
                                "type": "atom",
                                "value": 0.8,
                                "targets": [1, 2, 3]
                            }
                        ],
                        "circles": [
                            {
                                "type": "atom",
                                "targets": [1, 2],
                                "color": "#FF0000"
                            }
                        ],
                        "backbone_color": "#666666"
                    }
                }
            ]
        } 

def from_spec(spec: Union[dict, XenopictSpec]) -> List[str]:
    """Convert a specification into a list of SVG images.
    
    Args:
        spec: Either a dict matching the XenopictSpec schema or a XenopictSpec instance
        
    Returns:
        List of SVG strings, one for each molecule in the specification
    """
    if isinstance(spec, dict):
        spec = XenopictSpec(**spec)
    # TODO: Implement actual conversion to SVG
    return []

def from_json(json_str: str) -> List[str]:
    """Convert a JSON string specification into a list of SVG images.
    
    Args:
        json_str: JSON string matching the XenopictSpec schema
        
    Returns:
        List of SVG strings, one for each molecule in the specification
    """
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