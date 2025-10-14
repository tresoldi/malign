.PHONY: help quality format test test-cov test-fast bump-version build build-release clean install install-dev docs docs-clean bench

.DEFAULT_GOAL := help

PYTHON := python3
PIP := $(PYTHON) -m pip
TYPE ?= patch  # For version bumping

help: ## Show this help message
	@echo "MAlign Development Commands"
	@echo "==========================="
	@echo ""
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | \
	    awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

quality: ## Run code quality checks (ruff format --check, ruff check, mypy)
	ruff format --check .
	ruff check .
	mypy malign/ tests/ scripts/

format: ## Auto-format code with ruff
	ruff format .

test: ## Run test suite
	pytest tests/

test-cov: ## Run tests with coverage (HTML report)
	pytest --cov=malign --cov-report=html:tests/htmlcov \
	       --cov-report=term-missing tests/

test-fast: ## Run tests in parallel
	pytest -n auto tests/

bump-version: ## Bump version (TYPE=patch|minor|major), commit, and tag
	@CURRENT=$$(grep -o "__version__ = \"[^\"]*\"" malign/__init__.py | cut -d'"' -f2); \
	echo "==> Current version: $$CURRENT"; \
	IFS='.' read -r major minor patch <<< "$$CURRENT"; \
	if [ "$(TYPE)" = "major" ]; then NEW="$$((major + 1)).0.0"; \
	elif [ "$(TYPE)" = "minor" ]; then NEW="$$major.$$((minor + 1)).0"; \
	elif [ "$(TYPE)" = "patch" ]; then NEW="$$major.$$minor.$$((patch + 1))"; \
	else echo "Error: TYPE must be patch, minor, or major"; exit 1; fi; \
	echo "==> Bumping $(TYPE) version to $$NEW..."; \
	sed -i "s/__version__ = \"$$CURRENT\"/__version__ = \"$$NEW\"/" malign/__init__.py; \
	echo "⚠️  Please update CHANGELOG.md manually before committing!"; \
	read -p "Press Enter to commit and tag, or Ctrl+C to cancel..."; \
	git add malign/__init__.py; \
	git commit -m "chore: bump version to $$NEW"; \
	git tag -a "v$$NEW" -m "Release v$$NEW"; \
	echo "✓ Version bumped to $$NEW and tagged!"

build: ## Build package (creates dist/)
	$(PYTHON) -m build

build-release: clean quality test build ## Full release build (clean → quality → test → build)
	@echo "✓ Release build complete!"
	@ls -lh dist/

clean: ## Remove build artifacts, caches, and coverage reports
	rm -rf dist/ build/ *.egg-info
	rm -rf .coverage htmlcov/ tests/htmlcov/ coverage.xml
	rm -rf .pytest_cache .ruff_cache .mypy_cache
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true

install: ## Install package in development mode
	$(PIP) install -e .

install-dev: ## Install package with development dependencies
	$(PIP) install -e .[dev]

docs: ## Generate HTML documentation from Nhandu tutorial sources
	@for f in docs/tutorial_*.py; do \
		if [ -f "$$f" ]; then \
			echo "Generating $$(basename $$f .py).html..."; \
			nhandu "$$f" --format html -o "docs/$$(basename $$f .py).html"; \
		fi; \
	done

docs-clean: ## Remove generated HTML documentation
	rm -f docs/tutorial_*.html

bench: ## Run quick performance benchmarks
	$(PYTHON) scripts/benchmarks.py --quick --output benchmark_results_quick
