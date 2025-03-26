"""Layout specifications for xenopict's declarative API.

This module contains Pydantic models for specifying layout properties of molecule
visualizations, such as position and size.
"""

from pydantic import BaseModel, Field, ConfigDict


class Position(BaseModel):
    """Specifies a position in 2D space."""
    model_config = ConfigDict(frozen=True)

    x: float = Field(default=0.0, description="X coordinate")
    y: float = Field(default=0.0, description="Y coordinate")


class Size(BaseModel):
    """Specifies dimensions in 2D space."""
    model_config = ConfigDict(frozen=True)

    width: float = Field(default=100.0, description="Width in pixels")
    height: float = Field(default=100.0, description="Height in pixels")


class LayoutSpec(BaseModel):
    """Specifies layout properties for a molecule visualization."""
    model_config = ConfigDict(frozen=True)

    position: Position = Field(
        default_factory=Position,
        description="Position of the molecule visualization"
    )
    size: Size = Field(
        default_factory=Size,
        description="Size of the molecule visualization"
    )
    margin: float = Field(
        default=10.0,
        description="Margin around the molecule visualization in pixels"
    )
    padding: float = Field(
        default=5.0,
        description="Padding inside the molecule visualization in pixels"
    ) 