# XenoPict Design

XenoPict is a Python library for creating consistent and beautiful 2D molecular depictions. It provides both an imperative API (current) and a declarative API (in development) for molecular visualization tasks.

## Core Concepts

### Molecules and Depiction
- Uses RDKit as the underlying engine for molecular operations
- Focuses on 2D depiction and alignment
- Maintains consistent visual representation across multiple molecules
- Preserves chemical validity while optimizing visual layout

### Alignment System
The alignment system is built around three key concepts:
1. **Template Molecules**: Reference structures that define desired layouts
2. **Source Molecules**: Molecules to be aligned to templates
3. **Atom Pairs**: Explicit mappings between source and template atoms

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

## Planned Architecture (Declarative API)

### Overview
The declarative API will allow users to specify molecular layouts using a JSON-like syntax, inspired by VEGA's approach to data visualization.

### Schema Design
```json
{
  "version": "1.0",
  "molecules": [
    {
      "id": "mol1",
      "smiles": "CCO",
      "role": "template"
    },
    {
      "id": "mol2",
      "smiles": "CCCO",
      "alignTo": "mol1",
      "alignmentHints": {
        "type": "atomPairs",
        "pairs": [[2, 1]]  // align OH groups
      }
    }
  ],
  "style": {
    "bondLength": 1.5,
    "atomLabels": true,
    "colorScheme": "default"
  }
}
```

### Key Components (Planned)

#### 1. Schema Validation
- JSON Schema validation for structural correctness
- Chemical validation for molecular specifications
- Runtime type checking for Python interfaces

#### 2. Parser/Interpreter
- Converts declarative specifications into execution plan
- Resolves references and variables
- Handles template expansion
- Validates chemical semantics

#### 3. Execution Engine
- Builds molecule objects from specifications
- Applies alignments in optimal order
- Handles style application
- Manages coordinate systems

#### 4. Component System
- Reusable templates for common patterns
- Support for molecular fragments
- Composition of complex layouts
- Style inheritance and overrides

## Design Principles

### 1. Chemical Validity
- All operations preserve chemical validity
- Alignment never modifies molecular structure
- Clear separation between chemical and visual properties

### 2. Consistency
- Consistent visual representation across molecules
- Predictable alignment behavior
- Stable layouts for similar inputs

### 3. Flexibility
- Multiple approaches to alignment
- Support for both automatic and manual control
- Extensible style system

### 4. User Experience
- Sensible defaults for common cases
- Clear error messages
- Comprehensive documentation
- Progressive disclosure of complexity

## Implementation Notes

### Performance Considerations
- Lazy evaluation where possible
- Caching of intermediate results
- Efficient graph algorithms for multi-molecule alignment
- Smart defaults to minimize computation

### Error Handling
- Validation at schema level
- Chemical validity checks
- Clear error messages with suggestions
- Graceful fallbacks where appropriate

### Testing Strategy
1. Unit Tests
   - Individual component functionality
   - Edge cases and error conditions
   - Chemical validity preservation

2. Integration Tests
   - End-to-end workflows
   - Complex molecule sets
   - Style application

3. Performance Tests
   - Large molecule sets
   - Complex alignment scenarios
   - Memory usage

## Future Directions

### Planned Enhancements
1. 3D alignment support
2. Interactive visualization tools
3. Batch processing capabilities
4. Additional style options
5. Template library

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