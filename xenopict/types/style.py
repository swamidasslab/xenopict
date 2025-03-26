"""Style-related type definitions for xenopict's declarative API."""

from typing import Dict, List, Optional, Union, Literal
from pydantic import BaseModel, Field, validator
import re

class ColorScheme(BaseModel):
    """Color scheme configuration."""
    background: Optional[str] = Field(
        None,
        description="Background color (hex or name)"
    )
    atoms: Dict[str, str] = Field(
        default_factory=dict,
        description="Mapping of atom symbols to colors"
    )
    bonds: str = Field(
        "#000000",
        description="Default bond color"
    )
    highlights: str = Field(
        "#FF0000",
        description="Color for highlighted elements"
    )

    @validator('background', 'bonds', 'highlights', 'atoms')
    def validate_colors(cls, v):
        """Validate color format."""
        if isinstance(v, dict):
            for color in v.values():
                if not re.match(r'^#[0-9a-fA-F]{6}$', color):
                    raise ValueError(f"Invalid color format: {color}")
            return v
        if v is not None and not re.match(r'^#[0-9a-fA-F]{6}$', v):
            raise ValueError(f"Invalid color format: {v}")
        return v

class FontSpec(BaseModel):
    """Font configuration."""
    family: str = Field("Arial", description="Font family name")
    size: float = Field(12.0, description="Base font size")
    weight: Literal["normal", "bold"] = Field(
        "normal",
        description="Font weight"
    )

class StyleSpec(BaseModel):
    """Style configuration for molecule visualization."""
    colors: ColorScheme = Field(
        default_factory=ColorScheme,
        description="Color scheme configuration"
    )
    font: FontSpec = Field(
        default_factory=FontSpec,
        description="Font configuration"
    )
    bond_width: float = Field(
        1.0,
        description="Base bond width"
    )
    atom_size: float = Field(
        1.0,
        description="Base atom size multiplier"
    )
    show_hydrogens: bool = Field(
        False,
        description="Whether to show explicit hydrogens"
    )
    show_atom_numbers: bool = Field(
        False,
        description="Whether to show atom numbers"
    )
    show_electron_pairs: bool = Field(
        True,
        description="Whether to show electron pairs"
    )

    class Config:
        """Pydantic model configuration."""
        json_schema_extra = {
            "examples": [
                {
                    "colors": {
                        "background": "#FFFFFF",
                        "atoms": {"C": "#000000", "O": "#FF0000", "N": "#0000FF"},
                        "bonds": "#000000",
                        "highlights": "#FF0000"
                    },
                    "font": {
                        "family": "Arial",
                        "size": 12,
                        "weight": "normal"
                    },
                    "bond_width": 1.0,
                    "atom_size": 1.0,
                    "show_hydrogens": False,
                    "show_atom_numbers": True
                }
            ]
        } 