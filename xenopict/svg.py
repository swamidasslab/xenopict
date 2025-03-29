"""SVG type definition for xenopict."""

from typing import NewType

# SVG is just a string type, but we give it a distinct type for type checking
SVG = NewType("SVG", str)
