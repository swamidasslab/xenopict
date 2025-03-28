"""Type definitions for the declarative API.

This module contains the Pydantic models that define the schema for the declarative API.
The schema is used to validate input and generate JSON schema documentation.
"""

from typing import List, Union, Optional
from pydantic import BaseModel, Field, ConfigDict


class MoleculeSpec(BaseModel):
    """Specification for a single molecule.
    
    Examples:
        >>> # Basic molecule with just SMILES
        >>> spec = MoleculeSpec(smiles="CCO")
        >>> spec.smiles
        'CCO'
        
        >>> # Molecule with optional ID
        >>> spec = MoleculeSpec(smiles="CCO", id="ethanol")
        >>> spec.id
        'ethanol'
    """
    smiles: str = Field(
        ...,  # ... means required
        description="SMILES string representing the molecule",
        examples=["CCO", "c1ccccc1"],
    )
    id: Optional[str] = Field(
        None,
        description="Optional identifier for the molecule, used for reference in error messages",
        examples=["ethanol", "mol1"],
    )


class XenopictSpec(BaseModel):
    """Root specification for xenopict visualizations.
    
    This class represents the top-level schema for xenopict visualizations.
    It can contain either a single molecule or a list of molecules.
    
    Examples:
        >>> # Single molecule
        >>> spec = XenopictSpec(molecules=MoleculeSpec(smiles="CCO"))
        >>> isinstance(spec.molecules, MoleculeSpec)
        True
        
        >>> # Multiple molecules with alignment
        >>> spec = XenopictSpec(
        ...     molecules=[
        ...         MoleculeSpec(smiles="CCO"),
        ...         MoleculeSpec(smiles="CCCO")
        ...     ],
        ...     align=True
        ... )
        >>> isinstance(spec.molecules, list)
        True
        >>> len(spec.molecules)
        2
    """
    molecules: Union[MoleculeSpec, List[MoleculeSpec]] = Field(
        ...,
        description="One or more molecules to visualize",
        examples=[
            {"smiles": "CCO"},
            [{"smiles": "CCO"}, {"smiles": "CCCO"}]
        ],
    )
    align: Optional[bool] = Field(
        True,
        description="Whether to automatically align multiple molecules (default: True)",
        examples=[True, False],
    )

    model_config = ConfigDict(
        title="Xenopict Specification",
        json_schema_extra={
            "description": "Schema for defining molecule visualizations in xenopict",
            "examples": [
                {
                    "molecules": {"smiles": "CCO"},
                    "align": True
                },
                {
                    "molecules": [
                        {"smiles": "CCO", "id": "ethanol"},
                        {"smiles": "CCCO", "id": "propanol"}
                    ],
                    "align": True
                }
            ]
        }
    ) 