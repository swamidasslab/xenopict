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
	@echo "Downloading elkjs-svg..."
	curl -L https://cdn.jsdelivr.net/npm/elkjs-svg@latest/elkjs-svg.js -o $(JS_DIR)/elkjs-svg.js
	@if [ ! -s $(JS_DIR)/elkjs-svg.js ]; then \
		echo "Error: Failed to download elkjs-svg or file is empty"; \
		rm -f $(JS_DIR)/elkjs-svg.js; \
		exit 1; \
	fi
	@echo "Successfully downloaded elkjs-svg"
	@echo "Downloading elkjs-svg XML helpers..."
	mkdir -p $(JS_DIR)/helpers
	curl -L https://cdn.jsdelivr.net/npm/elkjs-svg@latest/helpers/xml.js -o $(JS_DIR)/helpers/xml.js
	@if [ ! -s $(JS_DIR)/helpers/xml.js ]; then \
		echo "Error: Failed to download XML helpers or file is empty"; \
		rm -f $(JS_DIR)/helpers/xml.js; \
		exit 1; \
	fi
	@echo "Successfully downloaded XML helpers"
	@echo "Processing elkjs-svg for browser compatibility..."
	python scripts/process_elkjs_svg.py
	@echo "Successfully processed elkjs-svg"

# Update all dependencies
deps: js-deps

# Clean up generated files
clean:
	rm -rf $(SCHEMA_DIR) $(COVERAGE_DIR) .coverage .pytest_cache __pycache__ build dist *.egg-info
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name "*.egg-info" -exec rm -r {} +
	find . -type f -name "*.pyc" -delete
