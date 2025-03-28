# XenoPict TODO List

## High Priority Tasks

### Immediate Fixes Needed
- [ ] Fix coordinate dimensionality issues
  - [ ] Ensure GetCoords consistently returns 2D coordinates
  - [ ] Update tests to expect 2D coordinates
  - [ ] Add validation to ensure 2D coordinates are generated
- [ ] Fix alignment consistency issues
  - [ ] Ensure consistent alignment between molecules
  - [ ] Fix hint-based alignment in auto_align_molecules
  - [ ] Add more robust alignment tests

### API Updates
- [x] Rename primary entry point from `create_molecules` to `parse`
- [x] Update all documentation to use `parse` as primary entry point
- [x] Add examples showing different ways to use `parse`
- [x] Add type hints and validation for all `parse` parameters
- [ ] Add migration guide for users moving from `create_molecules` to `parse`

### Declarative API Enhancements
- [x] Add support for molecule marking
  - [x] Individual atom marking with circles
  - [x] Substructure marking with bonds
  - [x] Proper SVG structure for marks
- [x] Add support for molecule styling
  - [x] Implement molecule colors
  - [x] Implement halo control
  - [x] Add tests for styling features
  - [x] Document styling options
- [ ] Add support for advanced features
  - [ ] Custom alignment methods
  - [ ] Style specifications
  - [ ] Grouping and composition
- [ ] Enhance schema
  - [ ] Add support for references and variables
  - [ ] Add support for templates and reusable components
  - [ ] Add support for custom alignment parameters
- [ ] Improve error handling
  - [ ] Add more descriptive error messages
  - [ ] Add validation for chemical validity
  - [ ] Add suggestions for common mistakes

### Code Improvements
- [x] Fix linter errors in test_declarative.py
  - [x] Update helper functions to match implementation
  - [x] Fix SVG structure verification
- [ ] Fix remaining linter errors:
  - "rdFMCS" unknown attribute of module "rdkit.Chem"
  - Type argument issues with FixtureFunction
- [ ] Add more type safety
  - [ ] Add runtime checks for RDKit molecule validity
  - [ ] Add better error messages for invalid SMILES
  - [ ] Add type guards for molecule conversions

## Medium Priority Tasks

### Documentation
- [ ] Add a user guide or tutorial in docs/ directory
- [ ] Add more complex examples showing real-world use cases
- [ ] Document best practices for molecule alignment
- [ ] Add visualization examples showing before/after alignments
- [ ] Add migration guide from imperative API to declarative API

### Testing
- [x] Add comprehensive tests for declarative API marking
  - [x] Individual atom marking
  - [x] Substructure marking
  - [x] Combined marking types
  - [x] SVG structure verification
- [x] Add comprehensive tests for molecule styling
  - [x] Color specification tests
  - [x] Halo control tests
  - [x] SVG style verification
  - [x] Error handling for invalid colors
- [ ] Add more comprehensive tests for declarative API
  - [ ] Edge cases for molecule specifications
  - [ ] Error handling tests
  - [ ] Integration tests with RDKit
  - [ ] Performance benchmarks
- [ ] Add stress tests for complex molecular structures
- [ ] Add more edge cases to test suite
- [ ] Add tests for coordinate dimensionality
- [ ] Add tests for alignment consistency

### Infrastructure
- [ ] Set up continuous integration
- [ ] Add code coverage reporting
- [ ] Add automated style checking
- [ ] Set up automated documentation building

## Lower Priority Tasks

### Code Improvements
- [ ] Consider adding type hints for networkx Graph objects
- [ ] Add more error checking for molecule validity
- [ ] Consider adding progress tracking for large molecule sets
- [ ] Add support for 3D alignment
- [ ] Add support for alignment quality metrics
- [ ] Consider adding visualization tools
- [ ] Add batch processing capabilities for large sets

## Notes
- [x] Initial version of declarative API is complete and working
- [x] Core functionality (molecule specification and alignment) is implemented
- [x] Primary entry point renamed to `parse` for clarity and consistency
- [x] Molecule marking functionality implemented and tested
- [ ] Focus should now be on fixing coordinate and alignment issues
- [ ] Consider gathering user feedback on current API before adding more features
- [ ] Maintain backward compatibility while adding new features
- [ ] Consider adding examples and documentation for common use cases
- [ ] Need to ensure consistent 2D coordinate handling throughout the codebase

## Completed Tasks

### Documentation
- [x] Added comprehensive module-level documentation for `alignment.py`
  - Emphasized `auto_align_molecules` as primary entrypoint
  - Added clear examples and use cases
  - Structured documentation to guide users to preferred methods
- [x] Documented `Alignment` class and all its methods
  - Added detailed docstrings with examples
  - Included validation and error cases
  - Added method chaining examples
- [x] Added doctests for all alignment methods
  - Basic functionality tests
  - Edge cases and error conditions
  - Coordinate verification tests

### Code Cleanup
- [x] Removed redundant `auto_alignment` function
- [x] Fixed networkx import and usage
- [x] Removed redundant doctests where high-level functions wrap lower-level ones
- [x] Removed "low-level" notes from core Alignment class methods
- [x] Moved types to declarative package
- [x] Removed unused type specifications
- [x] Cleaned up and simplified type hierarchy
- [x] Renamed primary entry point to `parse` for clarity and consistency
- [x] Removed deprecated `shaded_svg` function and all references to it
- [x] Fixed test helper functions to match implementation
  - [x] Updated _get_marked_bonds to look in mark group
  - [x] Fixed SVG structure verification in tests

### Declarative API Development
- [x] Design core schema for declarative molecule rendering
  - [x] Define JSON schema for molecular layout specifications
  - [x] Support for alignment specifications
  - [x] Simplified schema focused on core functionality
- [x] Implement schema validation
  - [x] Add JSON schema validation using Pydantic
  - [x] Add runtime type checking
  - [x] Add semantic validation for chemical structures
- [x] Create parser/interpreter for declarative specifications
  - [x] Implement basic parsing system
  - [x] Support for single and multiple molecules
  - [x] Support for automatic alignment
- [x] Add high-level convenience functions
  - [x] Primary entry point `parse` for all input types
  - [x] Support for both dict and model inputs
- [x] Documentation for declarative API
  - [x] Basic examples in docstrings
  - [x] Type hints and validation messages
  - [x] JSON schema documentation
- [x] Implement molecule marking functionality
  - [x] Individual atom marking with circles
  - [x] Substructure marking with bonds
  - [x] SVG structure for marks
  - [x] Comprehensive test coverage
- [x] Add molecule marking feature
- [x] Add hierarchical styling support (color and halo)
  - [x] Add style fields to XenopictSpec for defaults
  - [x] Allow MoleculeSpec to override defaults
  - [x] Implement style inheritance logic
  - [x] Add comprehensive tests for styling
- [x] Improve CI/CD
  - [x] Add schema generation command
  - [x] Add test commands
  - [x] Add coverage reporting
  - [x] Update documentation with development instructions

# In Progress

# TODO

- [ ] Add color validation to ensure valid color values
- [ ] Consider adding more style options (e.g., font size, stroke width)
- [ ] Add documentation and examples for styling in README.md
- [ ] Consider adding style presets/themes at XenopictSpec level
- [ ] Improve test coverage
  - [ ] Increase coverage in _version.py (51%)
  - [ ] Increase coverage in magic.py (55%)
  - [ ] Increase coverage in monkey.py (56%)
  - [ ] Increase coverage in drawer.py (71%)
- [ ] Add schema validation to ensure JSON input matches schema
- [ ] Consider adding schema versioning for future compatibility
