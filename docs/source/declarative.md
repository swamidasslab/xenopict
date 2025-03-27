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

## Molecule Alignment

The API supports several methods for aligning molecules:

### Using Atom Maps

You can align molecules by specifying atom mappings. There are two ways to specify mappings:

1. Using SMILES strings with atom mapping IDs:
```python
spec = {
    "molecule": [
        {
            "id": "ethanol",
            "smiles": "CCO",  # Structure without map IDs
            "map": "[CH3:1][CH2:2][OH:3]"  # Default mapping
        },
        {
            "id": "ethylamine",
            "smiles": "CCN",  # Structure without map IDs
            "map": "[CH3:1][CH2:2][NH2:3]",  # Default mapping
            "alignment": {
                "to": "ethanol",
                "method": "mapids"  # Uses default maps from both molecules
            }
        }
    ]
}
```

2. Using integer lists (more concise):
```python
spec = {
    "molecule": [
        {
            "id": "ethanol",
            "smiles": "CCO",
            "map": [1, 2, 3]  # Each index maps to target atom index
        },
        {
            "id": "ethylamine",
            "smiles": "CCN",
            "alignment": {
                "to": "ethanol",
                "method": "mapids",
                "map": [1, 2, 3]  # Maps atoms to corresponding indices
            }
        }
    ]
}
```

You can also specify partial mappings by using 0 for unmapped atoms:
```python
spec = {
    "molecule": [
        {
            "id": "ethanol",
            "smiles": "CCO",
            "map": [1, 2, 0]  # Third atom (O) is unmapped
        },
        {
            "id": "propanol",
            "smiles": "CCCO",
            "alignment": {
                "to": "ethanol",
                "method": "mapids",
                "map": [1, 2, 0, 0]  # Last two atoms unmapped
            }
        }
    ]
}
```

### Overriding Default Maps

Each molecule can have a default mapping (`map`), but you can override it for specific alignments:
```python
spec = {
    "molecule": [
        {
            "id": "mol1",
            "smiles": "CCO"
        },
        {
            "id": "mol2",
            "smiles": "CCN",
            "alignment": {
                "to": "mol1",
                "method": "mapids",
                "map": [1, 2, 3],  # Override mapping for this molecule
                "to_map": [1, 2, 3]  # Override mapping for target molecule
            }
        }
    ]
}
```

### Other Alignment Methods

#### Substructure Alignment
```python
spec = {
    "molecule": [
        {
            "id": "benzene",
            "smiles": "c1ccccc1"
        },
        {
            "id": "phenol",
            "smiles": "Oc1ccccc1",
            "alignment": {
                "to": "benzene",
                "method": "substructure",
                "params": {
                    "smarts": "c1ccccc1"  # Benzene ring pattern
                }
            }
        }
    ]
}
```

#### Manual Alignment
```python
spec = {
    "molecule": [
        {
            "id": "mol1",
            "smiles": "CCO"
        },
        {
            "id": "mol2",
            "smiles": "CCN",
            "alignment": {
                "to": "mol1",
                "method": "manual",
                "params": {
                    "atom_pairs": [[0, 0], [1, 1]]  # Align specific atom pairs
                }
            }
        }
    ]
}
```

#### Automatic Alignment
```python
spec = {
    "molecule": [
        {
            "id": "aspirin",
            "smiles": "CC(=O)Oc1ccccc1C(=O)O"
        },
        {
            "id": "salicylic",
            "smiles": "O=C(O)c1ccccc1O",
            "alignment": {
                "to": "aspirin",
                "method": "auto"  # Let RDKit determine best alignment
            }
        }
    ]
}
```

## Annotations

You can add various annotations to molecules:

### Shading
```python
spec = {
    "molecule": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "shading": [
            {
                "type": "atom",
                "value": 0.8,
                "targets": [1, 2, 3]
            },
            {
                "type": "bond",
                "value": 0.5,
                "targets": [[1, 2], [2, 3]]
            },
            {
                "type": "substructure",
                "value": 0.7,
                "targets": "c1ccccc1"  # SMARTS pattern
            }
        ]
    }
}
```

### Circles
```python
spec = {
    "molecule": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "circles": [
            {
                "type": "atom",
                "targets": [1, 2],
                "color": "#FF0000",
                "width": 2.0
            },
            {
                "type": "bond",
                "targets": [[3, 4]],
                "color": "#0000FF"
            }
        ]
    }
}
```

### Backbone Color
```python
spec = {
    "molecule": {
        "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
        "backbone_color": "#666666"
    }
}
``` 