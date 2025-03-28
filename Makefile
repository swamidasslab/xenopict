.PHONY: upload test coverage schema clean

# Directory structure
SCHEMA_DIR := docs/schema
COVERAGE_DIR := coverage_report

# Create necessary directories
$(SCHEMA_DIR):
	mkdir -p $(SCHEMA_DIR)

$(COVERAGE_DIR):
	mkdir -p $(COVERAGE_DIR)

upload:
	python setup.py build sdist
	twine upload dist/*
	git push --tags

# Generate JSON schema from Pydantic models
schema: $(SCHEMA_DIR)
	python scripts/generate_schema.py $(SCHEMA_DIR)

# Run tests
test:
	pytest xenopict/tests/ -v

# Run tests with coverage
coverage: $(COVERAGE_DIR)
	pytest xenopict/tests/ -v --cov=xenopict --cov-report=html:$(COVERAGE_DIR) --cov-report=term-missing

# Display coverage report in browser
show-coverage: coverage
	python -m webbrowser "$(COVERAGE_DIR)/index.html"

# Clean up generated files
clean:
	rm -rf $(SCHEMA_DIR) $(COVERAGE_DIR) .coverage .pytest_cache __pycache__ build dist *.egg-info
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name "*.egg-info" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
