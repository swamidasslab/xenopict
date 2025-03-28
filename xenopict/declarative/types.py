"""Type definitions for the declarative API.

This module contains the Pydantic models that define the schema for the declarative API.
The schema is used to validate input and generate JSON schema documentation.
"""

from typing import List, Optional, Union

from pydantic import BaseModel, ConfigDict, Field, model_validator


class MarkSpec(BaseModel):
    """Specification for marking parts of a molecule.

    This class defines how to highlight or mark specific parts of a molecule using
    Xenopict's marking capabilities. All indices are 0-based.

    There are two ways to mark a molecule:
    1. Mark individual atoms with circles using `atoms`
    2. Mark a substructure (defined by atoms and optional bonds) using `substructure_atoms`
       and optionally `substructure_bonds`

    JSON Examples:

    Mark specific atoms with circles:

    .. code-block:: json

        {
            "smiles": "CCO",
            "mark": {
                "atoms": [0, 1]  # Marks C and C atoms
            }
        }

    Mark a substructure (all connecting bonds included):

    .. code-block:: json

        {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [0, 1]  # Marks CC substructure
            }
        }

    Mark substructure with specific bonds:

    .. code-block:: json

        {
            "smiles": "CCO",
            "mark": {
                "substructure_atoms": [0, 1],
                "substructure_bonds": [[0, 1]]  # Only mark the C-C bond
            }
        }

    Python Examples:

    >>> # Create a MarkSpec for marking atoms
    >>> mark = MarkSpec(atoms=[0, 1])
    >>> mark.atoms
    [0, 1]
    >>> mark.substructure_atoms is None
    True

    >>> # Create a MarkSpec for marking a substructure
    >>> mark = MarkSpec(substructure_atoms=[0, 1], substructure_bonds=[(0, 1)])
    >>> mark.substructure_atoms
    [0, 1]
    >>> mark.substructure_bonds
    [(0, 1)]
    >>> mark.atoms is None
    True
    """

    atoms: Optional[list[int]] = Field(
        None, description="List of 0-based atom indices to mark with circles"
    )

    substructure_atoms: Optional[list[int]] = Field(
        None, description="List of 0-based atom indices defining a substructure to mark"
    )
    substructure_bonds: Optional[list[tuple[int, int]]] = Field(
        None,
        description="List of 0-based atom index pairs defining specific bonds to include in the substructure marking",
    )

    @model_validator(mode="after")
    def validate_mark_spec(model: "MarkSpec") -> "MarkSpec":
        """Validates that the marking specification is valid.

        Rules:
        1. If using substructure marking:
           - substructure_atoms must be provided if substructure_bonds is provided
           - substructure_bonds is optional
        2. At least one marking method must be specified (atoms or substructure)
        3. Cannot mix atom marking with substructure marking
        """
        # First validate substructure bonds require atoms
        if model.substructure_bonds is not None and model.substructure_atoms is None:
            raise ValueError(
                "Cannot specify 'substructure_bonds' without 'substructure_atoms'. "
                "To mark specific bonds, you must specify both the atoms and bonds "
                "that make up the substructure"
            )

        # Then validate bond indices if present
        if model.substructure_bonds is not None:
            if any(i < 0 or j < 0 for i, j in model.substructure_bonds):
                raise ValueError("Substructure bond indices must be non-negative")

            # Validate bond pairs are unique
            if len(model.substructure_bonds) != len(set(model.substructure_bonds)):
                raise ValueError("Duplicate bonds specified in substructure_bonds")

            # Validate bond pairs don't contain same atom
            if any(i == j for i, j in model.substructure_bonds):
                raise ValueError("Bond cannot connect an atom to itself")

        # Then validate atom indices if present
        if model.atoms is not None:
            if any(i < 0 for i in model.atoms):
                raise ValueError("Atom indices must be non-negative")

        if model.substructure_atoms is not None:
            if any(i < 0 for i in model.substructure_atoms):
                raise ValueError("Substructure atom indices must be non-negative")

        # Finally validate marking method requirements
        methods_used = sum(1 for x in [model.atoms, model.substructure_atoms] if x is not None)

        if methods_used == 0:
            raise ValueError(
                "Must specify at least one marking method: either 'atoms' for individual "
                "atom marking or 'substructure_atoms' for substructure marking"
            )

        if methods_used > 1:
            raise ValueError(
                "Cannot mix marking methods. Use either 'atoms' for individual atom marking "
                "or 'substructure_atoms' for substructure marking, not both"
            )

        return model


class MoleculeSpec(BaseModel):
    """Specification for a single molecule.

    Examples:
        >>> # Basic molecule with just SMILES
        >>> spec = MoleculeSpec(smiles="CCO")
        >>> spec.smiles
        'CCO'
        >>> spec.halo is None  # Defaults to None to use XenopictSpec value
        True

        >>> # Molecule with optional ID and explicitly disabled halo
        >>> spec = MoleculeSpec(smiles="CCO", id="ethanol", halo=False)
        >>> spec.id
        'ethanol'
        >>> spec.halo
        False

        >>> # Molecule with color override
        >>> spec = MoleculeSpec(smiles="CCO", color="blue")
        >>> spec.color
        'blue'
    """

    smiles: str = Field(
        ...,  # ... means required
        description="SMILES string representing the molecule",
        examples=["CCO", "c1ccccc1"],
    )
    mark: Optional[MarkSpec] = Field(
        None, description="Specification for marking parts of the molecule"
    )
    color: Optional[str] = Field(
        None,
        description="Color to use for the molecule's backbone and atom labels. If None, uses XenopictSpec color.",
        examples=["red", "blue", "green", "#000000"],
    )
    halo: Optional[bool] = Field(
        None,
        description="Whether to draw a subtle halo around the molecule. If None, uses XenopictSpec halo setting.",
        examples=[True, False],
    )
    id: Optional[str] = Field(
        None,
        description="Optional identifier for the molecule, used for reference in error messages",
        examples=["ethanol", "mol1"],
    )


class XenopictSpec(BaseModel):
    """Root specification for xenopict visualizations.

    This class represents the top-level schema for xenopict visualizations.
    It can contain either a single molecule or a list of molecules, and provides
    default styling options that apply to all molecules unless overridden.

    Examples:
        >>> # Single molecule with default styling
        >>> spec = XenopictSpec(
        ...     molecules=MoleculeSpec(smiles="CCO"),
        ...     color="red",  # All molecules will be red unless overridden
        ...     halo=True     # All molecules will have halos unless overridden
        ... )
        >>> isinstance(spec.molecules, MoleculeSpec)
        True
        >>> spec.color
        'red'
        >>> spec.halo
        True

        >>> # Multiple molecules with some overriding defaults
        >>> spec = XenopictSpec(
        ...     molecules=[
        ...         MoleculeSpec(smiles="CCO"),  # Uses default red color
        ...         MoleculeSpec(smiles="CCCO", color="blue")  # Overrides to blue
        ...     ],
        ...     color="red",
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
        examples=[{"smiles": "CCO"}, [{"smiles": "CCO"}, {"smiles": "CCCO"}]],
    )
    align: Optional[bool] = Field(
        True,
        description="Whether to automatically align multiple molecules (default: True)",
        examples=[True, False],
    )
    color: Optional[str] = Field(
        None,
        description="Default color to use for all molecules' backbones and atom labels",
        examples=["red", "blue", "green", "#000000"],
    )
    halo: bool = Field(
        True,
        description="Whether to draw subtle halos around molecules by default",
        examples=[True, False],
    )

    model_config = ConfigDict(
        title="Xenopict Specification",
        json_schema_extra={
            "description": "Schema for defining molecule visualizations in xenopict",
            "examples": [
                {"molecules": {"smiles": "CCO"}, "align": True, "color": "black", "halo": True},
                {
                    "molecules": [
                        {"smiles": "CCO", "id": "ethanol"},
                        {"smiles": "CCCO", "id": "propanol", "color": "blue"},
                    ],
                    "align": True,
                    "color": "red",
                    "halo": True,
                },
            ],
        },
    )
