# XenoPict TODO List

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

## High Priority Tasks

### Declarative API Development
- [ ] Design core schema for declarative molecule rendering
  - [ ] Define JSON schema for molecular layout specifications
  - [ ] Support for alignment specifications
  - [ ] Support for style specifications
  - [ ] Support for grouping and composition
- [ ] Implement schema validation
  - [ ] Add JSON schema validation
  - [ ] Add runtime type checking
  - [ ] Add semantic validation for chemical structures
- [ ] Create parser/interpreter for declarative specifications
  - [ ] Implement VEGA-like parsing system
  - [ ] Add support for references and variables
  - [ ] Add support for templates and reusable components
- [ ] Add high-level convenience functions
  - [ ] Factory functions for common use cases
  - [ ] Builder pattern for complex specifications
  - [ ] Fluent interface for chaining operations
- [ ] Documentation for declarative API
  - [ ] Tutorial with simple examples
  - [ ] Complex examples showing full capabilities
  - [ ] Best practices and patterns
  - [ ] Migration guide from imperative API

### Code Improvements
- [ ] Fix remaining linter errors:
  - "rdFMCS" unknown attribute of module "rdkit.Chem"
  - Type argument issues with FixtureFunction

## Medium Priority Tasks

### Documentation
- [ ] Add a user guide or tutorial in docs/ directory
- [ ] Add more complex examples showing real-world use cases
- [ ] Document best practices for molecule alignment
- [ ] Add visualization examples showing before/after alignments

### Testing
- [ ] Add comprehensive tests for declarative API
  - [ ] Schema validation tests
  - [ ] Parser tests
  - [ ] Integration tests with RDKit
  - [ ] Performance benchmarks
- [ ] Add stress tests for complex molecular structures
- [ ] Add more edge cases to test suite

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
- Primary focus should be on developing the declarative API while maintaining existing functionality
- Documentation should cover both imperative and declarative APIs
- Consider backward compatibility and migration path for existing users
- Gather user feedback on declarative API design before finalizing
