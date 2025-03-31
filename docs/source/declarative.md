---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.7
kernelspec:
  display_name: 'Python 3.9.13 (''python-3.9.13'': venv)'
  language: python
  name: python3
---

# Declarative API

The declarative API provides a simple way to depict molecules using a declarative schema. This is useful for applications that need to depict molecules without writing python code.

## Installation

The declarative API is included in the xenopict package.

```console
pip install xenopict
```

## Basic Usage

The declarative API is accessed through the `xenopict.declarative` module. The main function is `parse`, which takes a dictionary specification and returns a Xenopict object.

```{code-cell} ipython3
from xenopict.declarative import parse

parse({
    "molecules": {
        "smiles": "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"
    }
})
```

## Styling

The declarative API supports styling options for the molecule depiction.

```{code-cell} ipython3
parse({
    "molecules": {
        "smiles": "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
        "style": {
            "scale": 30,
            "add_atom_indices": True,
            "add_bond_indices": True
        }
    }
})
```

## Marking

The declarative API supports marking atoms and substructures.

```{code-cell} ipython3
parse({
    "molecules": {
        "smiles": "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl",
        "marks": [
            {
                "atoms": [0, 1, 2],
                "color": "red"
            },
            {
                "substructure_atoms": [7, 8, 9, 10],
                "color": "blue"
            }
        ]
    }
})
```

## Schema Validation

The declarative API validates all input against a schema. If the input is invalid, it raises a ValueError with a helpful message.

```{code-cell} ipython3
from xenopict.declarative import parse

try:
    parse({
        "molecules": {
            "smiles": "invalid smiles"
        }
    })
except ValueError as e:
    print(e)
```

Common validation errors include:

1. Invalid SMILES strings:

```{code-cell} ipython3
try:
    parse({
        "molecules": {
            "smiles": "not a valid smiles string"
        }
    })
except ValueError as e:
    print(e)
```

2. Invalid color values:

```{code-cell} ipython3
try:
    parse({
        "molecules": {
            "smiles": "CC",
            "marks": [
                {
                    "atoms": [0],
                    "color": "not a color"
                }
            ]
        }
    })
except ValueError as e:
    print(e)
```

3. Invalid marking specifications:

```{code-cell} ipython3
try:
    parse({
        "molecules": {
            "smiles": "CC",
            "marks": [
                {
                    "atoms": [2],  # atom index out of range
                    "color": "red"
                }
            ]
        }
    })
except ValueError as e:
    print(e)
```

## API Reference

For complete API details, see the [API Reference](api/xenopict.declarative.rst).

### Key Classes

- `XenopictSpec`: Root specification for visualizations
- `MoleculeSpec`: Specification for individual molecules
- `MarkSpec`: Specification for marking parts of molecules

### Functions

- `parse(spec)`: Convert a specification into a Xenopict object
- `parse_json(json_str)`: Parse a JSON string into a Xenopict object
- `parse_file(path)`: Parse a JSON file into a Xenopict object 