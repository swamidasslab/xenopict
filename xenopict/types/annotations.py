"""Annotation specifications for xenopict's declarative API.

This module contains Pydantic models for specifying annotations in molecule
visualizations, such as text, arrows, and highlights.
"""

from typing import List, Optional, Union
from pydantic import BaseModel, Field, ConfigDict

from .style import ColorScheme, FontSpec


class TextAnnotation(BaseModel):
    """Specifies a text annotation."""
    model_config = ConfigDict(frozen=True)

    text: str = Field(..., description="The text content of the annotation")
    x: float = Field(..., description="X coordinate of the text anchor point")
    y: float = Field(..., description="Y coordinate of the text anchor point")
    font: FontSpec = Field(
        default_factory=lambda: FontSpec(),
        description="Font properties for the text"
    )
    color: str = Field(
        default="black",
        description="Color of the text"
    )
    anchor: str = Field(
        default="middle",
        description="Text anchor position (start, middle, end)"
    )
    rotate: float = Field(
        default=0.0,
        description="Rotation angle in degrees"
    )


class ArrowAnnotation(BaseModel):
    """Specifies an arrow annotation."""
    model_config = ConfigDict(frozen=True)

    start_x: float = Field(..., description="X coordinate of the arrow start")
    start_y: float = Field(..., description="Y coordinate of the arrow start")
    end_x: float = Field(..., description="X coordinate of the arrow end")
    end_y: float = Field(..., description="Y coordinate of the arrow end")
    color: str = Field(
        default="black",
        description="Color of the arrow"
    )
    width: float = Field(
        default=1.0,
        description="Width of the arrow line"
    )
    head_size: float = Field(
        default=5.0,
        description="Size of the arrow head"
    )


class HighlightAnnotation(BaseModel):
    """Specifies a highlight annotation."""
    model_config = ConfigDict(frozen=True)

    x: float = Field(..., description="X coordinate of the highlight center")
    y: float = Field(..., description="Y coordinate of the highlight center")
    radius: float = Field(
        default=10.0,
        description="Radius of the highlight circle"
    )
    color: str = Field(
        default="yellow",
        description="Color of the highlight"
    )
    opacity: float = Field(
        default=0.3,
        description="Opacity of the highlight (0.0-1.0)"
    )


class AnnotationSpec(BaseModel):
    """Container for all annotations in a molecule visualization."""
    model_config = ConfigDict(frozen=True)

    texts: List[TextAnnotation] = Field(
        default_factory=list,
        description="List of text annotations"
    )
    arrows: List[ArrowAnnotation] = Field(
        default_factory=list,
        description="List of arrow annotations"
    )
    highlights: List[HighlightAnnotation] = Field(
        default_factory=list,
        description="List of highlight annotations"
    ) 