# Declarative API

Xenopict provides a declarative API for depicting molecules. This API offers a more functional and composable approach to molecular visualization using a schema-based specification.

## Installation

The declarative API is included in the main xenopict package:

```console
pip install xenopict
```

## Basic Usage

The declarative API uses a schema-based approach where you define your visualization using plain Python dictionaries and lists. Here's a simple example with a single molecule:

```python
from xenopict import parse, magic

# Create a specification for a single molecule
spec = {
    "molecules": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "id": "aspirin"
    }
}

# Convert spec into image(s)
parse(spec)  # Display the image
```

You can also specify multiple molecules:

```python
spec = {
    "molecules": [
        {
            "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
            "id": "aspirin"
        },
        {
            "smiles": "O=C(O)c1ccccc1O",
            "id": "salicylic_acid"
        }
    ]
}

parse(spec)
```

## Styling

The declarative API supports hierarchical styling where you can set default styles at the top level and override them for specific molecules.

### Default Styling

```python
spec = {
    "molecules": [
        {"smiles": "CCO"},
        {"smiles": "CCN"}
    ],
    "color": "blue",  # Default color for all molecules
    "halo": True     # Default halo setting for all molecules
}

parse(spec)
```

### Per-Molecule Styling

```python
spec = {
    "molecules": [
        {
            "smiles": "CCO",
            "color": "red"    # Override default color
        },
        {
            "smiles": "CCN",
            "halo": False    # Override default halo
        }
    ],
    "color": "blue",  # Default color
    "halo": True     # Default halo setting
}

parse(spec)
```

## Molecule Marking

You can mark specific parts of molecules using either individual atoms or substructures.

### Marking Individual Atoms

```python
spec = {
    "molecules": {
        "smiles": "CCO",
        "mark": {
            "atoms": [0, 1]  # Mark first two atoms
        }
    }
}

parse(spec)
```

### Marking Substructures

```python
spec = {
    "molecules": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "mark": {
            "substructure_atoms": [1, 2, 3],  # Mark atoms in substructure
            "substructure_bonds": [[1, 2], [2, 3]]  # Mark specific bonds
        }
    }
}

parse(spec)
```

## Schema Validation

The declarative API uses Pydantic models to validate input specifications. The schema is automatically generated and available in the documentation:

```python
# The schema is available in docs/schema/xenopict.json
```

### Common Validation Errors

1. Invalid SMILES string:
```python
spec = {
    "molecules": {
        "smiles": "invalid_smiles",  # Will raise ValidationError
    }
}

parse(spec)  # Will raise ValidationError
```

2. Invalid color value:
```python
spec = {
    "molecules": {
        "smiles": "CCO",
        "color": "not_a_color"  # Will raise ValidationError
    }
}

parse(spec)  # Will raise ValidationError
```

3. Invalid marking specification:
```python
spec = {
    "molecules": {
        "smiles": "CCO",
        "mark": {
            "atoms": [0, 1],
            "substructure_atoms": [1, 2]  # Can't mix marking methods
        }
    }
}

parse(spec)  # Will raise ValidationError
```

## API Reference

For complete API details, see the [API Reference](api/xenopict.declarative.html).

### Key Classes

- `XenopictSpec`: Root specification for visualizations
- `MoleculeSpec`: Specification for individual molecules
- `MarkSpec`: Specification for marking parts of molecules

### Functions

- `parse(spec)`: Convert a specification into a Xenopict object
- `parse_json(json_str)`: Parse a JSON string into a Xenopict object
- `parse_file(path)`: Parse a JSON file into a Xenopict object 