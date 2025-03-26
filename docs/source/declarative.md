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
from xenopict.declarative import from_spec

# Create a specification for a single molecule
spec = {
    "molecule": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
    }
}

# Convert spec into image(s)
images = from_spec(spec)  # Returns a list of SVG images
images[0]  # Display first image
```

You can also load specifications directly from JSON files or strings:

```python
from xenopict.declarative import from_json, from_file

# Load from JSON string
json_str = '''
{
    "molecule": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"
    }
}
'''
images = from_json(json_str)

# Load from JSON file
images = from_file('my_visualization.json')
```

Here's an example of a more complex visualization with annotations:

```python
spec = {
    "molecule": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "shading": [
            {
                "type": "atom",
                "value": 0.8,
                "targets": [1, 2, 3]
            }
        ],
        "circles": [
            {
                "type": "atom",
                "targets": [1, 2],
                "color": "#FF0000"
            }
        ],
        "backbone_color": "#666666"
    }
}

images = from_spec(spec)
```

The specification can include:
- A single molecule or list of molecules
- Atom and bond shading
- Circled atoms and bonds
- Molecule alignments
- Custom colors and styles

When working with multiple molecules that need to be aligned, you'll need to provide IDs for the molecules being referenced:

```python
spec = {
    "molecule": [
        {
            "id": "mol1",  # ID required for alignment reference
            "smiles": "[CH3:1][CH2:2][OH:3]"
        },
        {
            "id": "mol2",  # ID required for alignment reference
            "smiles": "[CH3:1][CH2:2][NH2:3]"
        }
    ],
    "alignments": [
        {
            "method": "mapids",
            "reference_mol": "mol1",
            "target_mol": "mol2"
        }
    ]
}

images = from_spec(spec)  # Returns list of two aligned molecule images
```

Note that molecule IDs are optional and only required when the molecule needs to be referenced elsewhere in the specification (such as in alignments). The `from_spec`, `from_json`, and `from_file` functions always return a list of SVG images, even when the input specification contains only a single molecule. 