.PHONY: upload test coverage schema clean deps

# Directory structure
SCHEMA_DIR := docs/schema
COVERAGE_DIR := coverage_report
JS_DIR := xenopict/layout/js

# Create necessary directories
$(SCHEMA_DIR):
	mkdir -p $(SCHEMA_DIR)

$(COVERAGE_DIR):
	mkdir -p $(COVERAGE_DIR)

$(JS_DIR):
	mkdir -p $(JS_DIR)

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
	COVERAGE_FILE=.coverage pytest xenopict/tests/ --cov=xenopict --cov-config=coverage.ini --cov-report=html:$(COVERAGE_DIR) --cov-report=xml --cov-report=json --cov-report=term-missing

# Display coverage report in browser
show-coverage: coverage
	python -m webbrowser "$(COVERAGE_DIR)/index.html"

# Download and install JavaScript dependencies
js-deps: $(JS_DIR)
	@echo "Downloading elkjs..."
	curl -L https://cdn.jsdelivr.net/npm/elkjs@0.10.0/lib/elk.bundled.js -o $(JS_DIR)/elk.js
	@if [ ! -s $(JS_DIR)/elk.js ]; then \
		echo "Error: Failed to download elkjs or file is empty"; \
		rm -f $(JS_DIR)/elk.js; \
		exit 1; \
	fi
	@echo "Successfully downloaded elkjs"

# Update all dependencies
deps: js-deps

# Clean up generated files
clean:
	rm -rf $(SCHEMA_DIR) $(COVERAGE_DIR) .coverage .pytest_cache __pycache__ build dist *.egg-info
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name "*.egg-info" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
