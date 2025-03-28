# XenoPict Design

XenoPict is a Python library for creating consistent and beautiful 2D molecular depictions. It provides both an imperative API (current) and a declarative API (in development) for molecular visualization tasks.

## Core Concepts

### Molecules and Depiction
- Uses RDKit as the underlying engine for molecular operations
- Focuses strictly on 2D depiction and alignment
- Maintains consistent visual representation across multiple molecules
- Preserves chemical validity while optimizing visual layout
- Enforces 2D coordinate handling throughout the codebase

### API Design Principles
- Clear, action-oriented function names that reflect their purpose
- Consistent with Python's "explicit is better than implicit"
- Follows established patterns from popular declarative APIs
- Primary entry point named `parse` because:
  - Explicit about transforming specifications into objects
  - Common in declarative APIs (e.g., Pydantic)
  - Implies both validation and transformation
  - Captures the full process: reading, validating, and creating objects
  - Aligns with Python's emphasis on readability and clarity

### Alignment System
The alignment system is built around three key concepts:
1. **Template Molecules**: Reference structures that define desired layouts
2. **Source Molecules**: Molecules to be aligned to templates
3. **Atom Pairs**: Explicit mappings between source and template atoms

### Coordinate System
- Strictly 2D coordinates for all molecular depictions
- Consistent coordinate handling across all operations
- Clear separation between 2D and 3D functionality
- Validation to ensure 2D coordinates are maintained

### Marking System
The marking system provides two ways to highlight parts of molecules:
1. **Atom Marking**: Individual atoms are marked with circles
2. **Substructure Marking**: Connected atoms and bonds are marked as a unit

Key features:
- Clear visual distinction between marking types
- Consistent SVG structure for all marks
- Support for multiple marks per molecule
- Proper layering in SVG output

## Current Architecture (Imperative API)

### Alignment Class
The core `Alignment` class represents a 2D alignment between two molecules:
```python
class Alignment:
    atom_pairs: list[tuple[int, int]]  # (source_idx, template_idx) pairs
    score: float                        # Alignment quality score
    source_mol: Chem.Mol               # Molecule being aligned
    template_mol: Chem.Mol             # Template to align to
```

### Factory Methods
Multiple ways to create alignments, catering to different use cases:
1. `from_mcs()`: Automatic alignment using maximum common substructure
2. `from_aligned_atoms()`: Manual alignment with explicit atom pairs
3. `from_mapids()`: Alignment using atom map IDs (useful for reactions)
4. `from_indices()`: Alignment using atom index mappings

### Multi-Molecule Alignment
The `auto_align_molecules()` function provides intelligent multi-molecule alignment:
- Builds a maximum spanning tree of alignments
- Uses molecule size and similarity to determine optimal alignment order
- Supports hint-based alignment for user control
- Ensures global consistency across all molecules
- Maintains 2D coordinate consistency throughout alignment process

## Declarative API

### Overview
The declarative API allows users to specify molecular layouts and markings using a JSON-like syntax, inspired by VEGA's approach to data visualization.

### Schema Design
The declarative API uses a simple, focused schema that prioritizes the core functionality of molecule visualization:

```json
{
  "align": true,              // Optional: Whether to align molecules (defaults to true)
  "molecules": [
    {
      "smiles": "CCO",       // Required: SMILES string of the molecule
      "mark": {              // Optional: Marking specification
        "atoms": [0, 1]      // Mark atoms 0 and 1 with circles
      },
      "color": "#FF0000",    // Optional: Color for molecule (default: black)
      "halo": true          // Optional: Draw halo around molecule (default: true)
    },
    {
      "smiles": "CCCO",      // Required: SMILES string of the molecule
      "mark": {              // Optional: Marking specification
        "substructure_atoms": [0, 1, 2],  // Mark atoms 0, 1, 2 as substructure
        "substructure_bonds": [[0, 1]]    // Only mark bond between atoms 0 and 1
      },
      "color": "blue",       // Colors can be names or hex codes
      "halo": false         // Disable halo for this molecule
    }
  ]
}
```

Key Features:
1. **Minimal Required Fields**: Only `smiles` is required for each molecule
2. **Simple Alignment**: Optional top-level `align` flag controls alignment of all molecules
3. **Automatic Alignment**: When enabled, uses maximum common substructure (MCS) by default
4. **Flexible Marking**: Support for both atom and substructure marking
5. **Clear SVG Structure**: Consistent organization of SVG elements
6. **Visual Styling**: Support for molecule colors and visibility-enhancing halos

Example Usage:
```python
from xenopict import parse

# Create molecules with different styles
spec = {
    "molecules": [
        {
            "smiles": "CCO",      # ethanol
            "color": "red",       # Use named color
            "mark": {
                "atoms": [0, 1]   # Mark C atoms with circles
            }
        },
        {
            "smiles": "CCCO",     # propanol
            "color": "#0000FF",   # Use hex color
            "halo": false,        # Disable halo
            "mark": {
                "substructure_atoms": [0, 1],  # Mark CC substructure
                "substructure_bonds": [[0, 1]]  # Only mark C-C bond
            }
        }
    ]
}

# Returns list of Xenopict objects, molecules automatically aligned
xenopicts = parse(spec)

# Convert to SVG strings if needed
svgs = [str(x) for x in xenopicts]
```

This schema:
- Makes the API more approachable
- Reduces cognitive load
- Focuses on common use cases
- Provides a solid foundation for future extensions

### Key Components

#### 1. Schema Validation
- JSON Schema validation using Pydantic
- Chemical validation for molecular specifications
- Runtime type checking for Python interfaces
- Validation for marking specifications

#### 2. Parser/Interpreter
- Converts declarative specifications into execution plan
- Resolves references and variables
- Handles template expansion
- Validates chemical semantics
- Processes marking specifications

#### 3. Execution Engine
- Builds molecule objects from specifications
- Applies alignments in optimal order
- Handles style application
- Manages coordinate systems
- Creates SVG marks in proper layers

#### 4. SVG Structure
The SVG output follows a consistent structure:
```xml
<svg>
  <g class="shading">...</g>
  <g class="mol_halo">
    <!-- Optional halo for better visibility -->
    <use href="#mol_halo_xeno_id" class="halo" .../>
  </g>
  <g class="lines" style="stroke:#000000">  <!-- Color can be customized -->
    <g id="lines_xeno_id">
      <path class="bond-0 atom-0 atom-1" .../>
    </g>
  </g>
  <g class="text">...</g>
  <g class="overlay">
    <g class="mark">...</g>
  </g>
</svg>
```

## Design Principles

### 1. Chemical Validity
- All operations preserve chemical validity
- Alignment never modifies molecular structure
- Clear separation between chemical and visual properties

### 2. Consistency
- Consistent visual representation across molecules
- Predictable alignment behavior
- Stable layouts for similar inputs
- Strict adherence to 2D coordinate system
- Consistent coordinate handling across all operations
- Consistent SVG structure and class naming

### 3. Flexibility
- Multiple approaches to alignment
- Support for both automatic and manual control
- Multiple marking types for different needs
- Clear coordinate system boundaries
- Extensible SVG structure

### 4. User Experience
- Sensible defaults for common cases
  - Black color for molecules by default
  - Halos enabled by default for better visibility
- Clear error messages
- Comprehensive documentation
- Progressive disclosure of complexity
- Predictable coordinate handling
- Intuitive marking system
- Flexible styling options

## Implementation Notes

### Performance Considerations
- Lazy evaluation where possible
- Caching of intermediate results
- Efficient graph algorithms for multi-molecule alignment
- Smart defaults to minimize computation
- Optimized 2D coordinate generation and manipulation

### Error Handling
- Validation at schema level
- Chemical validity checks
- Clear error messages for marking issues
- Helpful suggestions for common mistakes

## Testing Strategy
1. Unit Tests
   - Individual component functionality
   - Edge cases and error conditions
   - Chemical validity preservation
   - Coordinate system validation
   - Alignment consistency verification

2. Integration Tests
   - End-to-end workflows
   - Complex molecule sets
   - Style application
   - Multi-molecule alignment scenarios
   - Coordinate system consistency

3. Performance Tests
   - Large molecule sets
   - Complex alignment scenarios
   - Memory usage
   - Coordinate generation and manipulation efficiency

## Future Directions

### Planned Enhancements
1. Enhanced alignment algorithms
2. Interactive visualization tools
3. Improved coordinate system validation
4. Advanced alignment quality metrics
5. Better alignment hint handling

### Areas for Research
1. Machine learning for alignment hints
2. Novel layout algorithms
3. Automated style selection
4. Performance optimizations

## Migration Path
For users moving from imperative to declarative API:
1. Direct equivalents for common operations
2. Helper functions for conversion
3. Mixed-mode operation during transition
4. Comprehensive migration guide

## Styling

### Hierarchical Styling

Styling in xenopict follows a hierarchical model:

1. Default styles can be set at the `XenopictSpec` level:
   ```json
   {
     "molecules": [...],
     "color": "red",     // Default color for all molecules
     "halo": true       // Default halo setting for all molecules
   }
   ```

2. Individual molecules can override these defaults in their `MoleculeSpec`:
   ```json
   {
     "molecules": [
       {"smiles": "CCO"},                    // Uses default styles
       {"smiles": "CCCO", "color": "blue"},  // Overrides color only
       {"smiles": "CCCCO", "halo": false}    // Overrides halo only
     ],
     "color": "red",
     "halo": true
   }
   ```

3. Style inheritance:
   - If a molecule doesn't specify a style, it uses the `XenopictSpec` default
   - If `XenopictSpec` doesn't specify a color, no color is applied
   - If `XenopictSpec` doesn't specify a halo setting, halos are enabled by default

## Development

### JSON Schema

The declarative API's schema is automatically generated from the Pydantic models and stored in `docs/schema/xenopict.json`. This serves as the authoritative reference for the JSON format.

To generate the schema:
```bash
make schema
```

### Testing and Coverage

Tests are written using pytest and can be run with:
```bash
make test        # Run tests
make coverage    # Run tests with coverage report
```

Coverage reports are generated in HTML format in the `coverage_report/` directory and can be viewed with:
```bash
make show-coverage
``` 