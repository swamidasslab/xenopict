"""Molecule-related type definitions for xenopict's declarative API."""

from typing import List, Optional, Union
from pydantic import BaseModel, Field, validator
import re

class AtomSpec(BaseModel):
    """Specification for atom-level customization."""
    index: int = Field(..., description="0-based atom index in the molecule")
    label: Optional[str] = Field(None, description="Custom atom label")
    color: Optional[str] = Field(None, description="Custom atom color (hex or name)")
    size: Optional[float] = Field(None, description="Custom atom size multiplier")

    @validator('color')
    def validate_color(cls, v):
        """Validate color format."""
        if v is None:
            return v
        if re.match(r'^#[0-9a-fA-F]{6}$', v):
            return v
        # Add validation for named colors if needed
        return v

class BondSpec(BaseModel):
    """Specification for bond-level customization."""
    atoms: List[int] = Field(..., description="List of 2 atom indices defining the bond")
    style: Optional[str] = Field(None, description="Bond style (solid, dashed, etc)")
    width: Optional[float] = Field(None, description="Custom bond width multiplier")
    color: Optional[str] = Field(None, description="Custom bond color")

    @validator('atoms')
    def validate_atoms(cls, v):
        """Validate that exactly 2 atoms are specified."""
        if len(v) != 2:
            raise ValueError("Bond must be defined by exactly 2 atom indices")
        return v

class MoleculeSpec(BaseModel):
    """Specification for a single molecule visualization."""
    smiles: str = Field(..., description="SMILES string of the molecule")
    name: Optional[str] = Field(None, description="Name or identifier for the molecule")
    atoms: List[AtomSpec] = Field(
        default_factory=list,
        description="Custom atom specifications"
    )
    bonds: List[BondSpec] = Field(
        default_factory=list,
        description="Custom bond specifications"
    )
    style: Optional["StyleSpec"] = Field(
        None,
        description="Molecule-specific style overrides"
    )
    highlight_atoms: Optional[List[int]] = Field(
        None,
        description="List of atom indices to highlight"
    )
    highlight_bonds: Optional[List[List[int]]] = Field(
        None,
        description="List of atom index pairs defining bonds to highlight"
    )

    class Config:
        """Pydantic model configuration."""
        json_schema_extra = {
            "examples": [
                {
                    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
                    "name": "Aspirin",
                    "highlight_atoms": [1, 2, 3],
                    "atoms": [
                        {"index": 0, "label": "Me", "color": "#FF0000"}
                    ],
                    "bonds": [
                        {"atoms": [1, 2], "style": "dashed"}
                    ]
                }
            ]
        } 