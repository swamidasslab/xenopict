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
- [ ] Update all documentation to use `parse` as primary entry point
- [ ] Add migration guide for users moving from `create_molecules` to `parse`
- [ ] Add examples showing different ways to use `parse`
- [ ] Add type hints and validation for all `parse` parameters

### Declarative API Enhancements
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

### API Updates
- [ ] Update all documentation to use `parse` as primary entry point
- [ ] Add migration guide for users moving from `create_molecules` to `parse`
- [ ] Add examples showing different ways to use `parse`
- [ ] Add type hints and validation for all `parse` parameters
