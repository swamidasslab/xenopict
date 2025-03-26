"""Base type definitions for xenopict's declarative API."""

from typing import Dict, List, Optional, Union, Literal
from pydantic import BaseModel, Field, validator, confloat
from enum import Enum

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
    id: str = Field(..., description="Unique identifier for the molecule")
    smiles: str = Field(..., description="SMILES string of the molecule")
    shading: Optional[List[ShadeSpec]] = None
    circles: Optional[List[CircleSpec]] = None
    backbone_color: Optional[str] = Field(
        None, 
        pattern="^#[0-9a-fA-F]{6}$",
        description="Color for the molecule's backbone"
    )

class XenopictSpec(BaseModel):
    """Root specification for xenopict visualizations."""
    molecules: List[MoleculeSpec] = Field(
        ...,
        description="List of molecules to render"
    )
    alignments: Optional[List[AlignmentSpec]] = None

    @validator("molecules")
    def validate_molecules(cls, v):
        """Validate that there is at least one molecule."""
        if not v:
            raise ValueError("At least one molecule must be specified")
        return v

    class Config:
        """Pydantic model configuration with examples."""
        json_schema_extra = {
            "examples": [
                # Example 1: Simple multi-molecule display
                {
                    "molecules": [
                        {
                            "id": "aspirin",
                            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
                        },
                        {
                            "id": "ibuprofen",
                            "smiles": "CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O"
                        }
                    ]
                },
                
                # Example 2: Molecule with shading and circles
                {
                    "molecules": [
                        {
                            "id": "aspirin_annotated",
                            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                            "shading": [
                                {
                                    "type": "atom",
                                    "value": 0.8,
                                    "targets": [1, 2, 3]
                                },
                                {
                                    "type": "substructure",
                                    "value": -0.5,
                                    "targets": "c1ccccc1",
                                    "zero_centered": False
                                }
                            ],
                            "circles": [
                                {
                                    "type": "atom",
                                    "targets": [1, 2],
                                    "color": "#FF0000"
                                },
                                {
                                    "type": "bond",
                                    "targets": [[3, 4]],
                                    "color": "#0000FF",
                                    "width": 2.0
                                }
                            ],
                            "backbone_color": "#666666"
                        }
                    ]
                },
                
                # Example 3: Molecule alignment with SMARTS
                {
                    "molecules": [
                        {
                            "id": "ref_mol",
                            "smiles": "c1ccccc1NC(=O)O"
                        },
                        {
                            "id": "target_mol",
                            "smiles": "c1ccccc1NC(=O)CC"
                        }
                    ],
                    "alignments": [
                        {
                            "method": "substructure",
                            "reference_mol": "ref_mol",
                            "target_mol": "target_mol",
                            "params": {
                                "smarts": "c1ccccc1N"
                            }
                        }
                    ]
                },
                
                # Example 4: Map ID-based alignment
                {
                    "molecules": [
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
                }
            ]
        } 