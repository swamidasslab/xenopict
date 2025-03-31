"""XenoPict: A library for molecular visualization with a focus on 2D depiction."""

# ruff: noqa: I001

from ._version import __version__
from .drawer import Xenopict
from .declarative import parse
from . import magic  # noqa: F401

__all__ = ["Xenopict", "parse", "__version__"]
