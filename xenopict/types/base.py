"""Base type definitions for xenopict's declarative API."""

from typing import Dict, List, Optional, Union, Literal
from pydantic import BaseModel, Field, field_validator, confloat, ConfigDict
from enum import Enum
import json

class AlignmentMethod(str, Enum):
    """Methods for aligning molecules."""
    AUTO = "auto"  # Let RDKit automatically align
    SUBSTRUCTURE = "substructure"  # Align based on common substructure
    MAPIDS = "mapids"  # Align based on atom map IDs in SMILES
    MANUAL = "manual"  # Manual alignment using atom indices

# Type for map specifications - either SMILES with map IDs or list of integers
MapSpec = Union[str, List[int]]

class AlignmentSpec(BaseModel):
    """Specification for how a molecule should be aligned to another molecule."""
    to: str = Field(..., description="ID of the molecule to align to")
    method: AlignmentMethod = Field(
        default=AlignmentMethod.AUTO,
        description="Method to use for alignment"
    )
    map: Optional[MapSpec] = Field(
        None,
        description="Mapping for this molecule in this alignment. Can be either:\n"
                   "- SMILES string with atom mapping IDs (e.g. '[CH3:1][CH2:2][OH:3]')\n"
                   "- List of integers where each index maps to target atom index, 0 means unmapped\n"
                   "If not provided, uses molecule's default map"
    )
    to_map: Optional[MapSpec] = Field(
        None,
        description="Mapping for the target molecule in this alignment. Same format as 'map'.\n"
                   "If not provided, uses target molecule's default map"
    )
    params: Optional[Dict] = Field(
        None,
        description="Method-specific parameters",
        examples=[
            # For SUBSTRUCTURE
            {"smarts": "c1ccccc1"},  
            # For MANUAL
            {"atom_pairs": [[1, 1], [2, 2], [3, 3]]}  # List of [this_atom_idx, other_atom_idx]
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
    id: Optional[str] = Field(None, description="Optional identifier for the molecule, required only when referenced by alignments")
    smiles: Optional[str] = Field(None, description="SMILES string of the molecule structure (without atom mapping IDs)")
    smarts: Optional[str] = Field(None, description="SMARTS pattern of the molecule")
    map: Optional[MapSpec] = Field(
        None,
        description="Default mapping for this molecule, used for alignments if not overridden. Can be either:\n"
                   "- SMILES string with atom mapping IDs (e.g. '[CH3:1][CH2:2][OH:3]')\n"
                   "- List of integers where each index maps to target atom index, 0 means unmapped"
    )
    shading: Optional[List[ShadeSpec]] = None
    circles: Optional[List[CircleSpec]] = None
    backbone_color: Optional[str] = Field(
        None, 
        pattern="^#[0-9a-fA-F]{6}$",
        description="Color for the molecule's backbone"
    )
    alignment: Optional[AlignmentSpec] = Field(
        None,
        description="Specification for how this molecule should be aligned to another molecule"
    )

    @field_validator("map")
    @classmethod
    def validate_map(cls, v: Optional[MapSpec], info) -> Optional[MapSpec]:
        """Validate map format."""
        if v is not None:
            if isinstance(v, list):
                # Validate list of integers
                if not all(isinstance(x, int) and x >= 0 for x in v):
                    raise ValueError("Map list must contain non-negative integers")
            elif isinstance(v, str):
                # Validate SMILES with map IDs
                if not any(f":{i}" in v for i in range(1, len(v))):
                    raise ValueError("Map SMILES must contain atom mapping IDs (e.g. '[CH3:1]')")
            else:
                raise ValueError("Map must be either a SMILES string with atom mapping IDs or a list of integers")
        return v

    @field_validator("map")
    @classmethod
    def validate_map_smiles(cls, v: Optional[MapSpec], info) -> Optional[MapSpec]:
        """Validate that map is provided if molecule has alignments using mapids method."""
        if v is None and info.data.get("alignment") is not None:
            alignment = info.data["alignment"]
            if alignment.method == AlignmentMethod.MAPIDS and alignment.map is None:
                raise ValueError("Must provide either molecule.map or alignment.map when using mapids alignment method")
        return v

    @field_validator("alignment")
    @classmethod
    def validate_alignment_id(cls, v: Optional[AlignmentSpec], info) -> Optional[AlignmentSpec]:
        """Validate that molecule has ID if alignment is specified."""
        if v is not None and info.data.get("id") is None:
            raise ValueError("Molecule must have an ID if alignment is specified")
        return v

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

    @field_validator("molecule")
    @classmethod
    def validate_alignments(cls, v: Union[MoleculeSpec, List[MoleculeSpec]]) -> Union[MoleculeSpec, List[MoleculeSpec]]:
        """Validate that alignment references are valid."""
        molecules = v if isinstance(v, list) else [v]
        
        # Create a set of available IDs
        molecule_ids = {m.id for m in molecules if m.id is not None}
        
        # Check each molecule's alignment
        for mol in molecules:
            if mol.alignment is not None:
                if mol.alignment.to not in molecule_ids:
                    raise ValueError(f"Alignment target molecule '{mol.alignment.to}' not found")
                
                # Check for circular alignments
                seen = {mol.id}
                current = mol
                while current.alignment is not None:
                    target_id = current.alignment.to
                    if target_id in seen:
                        raise ValueError(f"Circular alignment detected involving molecule '{mol.id}'")
                    seen.add(target_id)
                    # Find the target molecule
                    current = next(m for m in molecules if m.id == target_id)
        
        return v

    model_config = ConfigDict(
        json_schema_extra={
            "examples": [
                # Example 1: Single molecule
                {
                    "molecule": {
                        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
                    }
                },
                
                # Example 2: Multiple molecules with alignment using default maps (SMILES style)
                {
                    "molecule": [
                        {
                            "id": "mol1",
                            "smiles": "CCO",  # Structure without map IDs
                            "map": "[CH3:1][CH2:2][OH:3]"  # Default mapping using SMILES
                        },
                        {
                            "id": "mol2",
                            "smiles": "CCN",  # Structure without map IDs
                            "map": "[CH3:1][CH2:2][NH2:3]",  # Default mapping using SMILES
                            "alignment": {
                                "to": "mol1",
                                "method": "mapids"
                            }
                        }
                    ]
                },
                
                # Example 3: Multiple molecules with alignment using specific maps (integer list style)
                {
                    "molecule": [
                        {
                            "id": "mol1",
                            "smiles": "CCO",  # Structure without map IDs
                            "map": [1, 2, 3]  # Default mapping using integers
                        },
                        {
                            "id": "mol2",
                            "smiles": "CCN",  # Structure without map IDs
                            "alignment": {
                                "to": "mol1",
                                "method": "mapids",
                                "map": [1, 2, 3],  # Specific mapping using integers
                                "to_map": [1, 2, 3]  # Target mapping using integers
                            }
                        }
                    ]
                },

                # Example 4: Partial mapping with unmapped atoms
                {
                    "molecule": [
                        {
                            "id": "mol1",
                            "smiles": "CCO",
                            "map": [1, 2, 0]  # Third atom (O) is unmapped
                        },
                        {
                            "id": "mol2",
                            "smiles": "CCCO",
                            "alignment": {
                                "to": "mol1",
                                "method": "mapids",
                                "map": [1, 2, 0, 0]  # Last two atoms unmapped
                            }
                        }
                    ]
                },
                
                # Example 5: Single annotated molecule
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
    )

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