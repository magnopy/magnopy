.DEFAULT_GOAL := help

# Sphinx-related variables
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = _build

help:
	@echo "\x1b[31m"
	@echo "Please specify what do you want to do."
	@echo "\x1b[0m"
	@echo "Available options are:\n"
	@echo "    help - show this message"
	@echo "    html - build the html docs"
	@echo "    clean-html - clean all files from docs and build html docs from scratch"
	@echo "    html-from-zero - as clean-html + replot all figures"
	@echo "    doctest - run doctests"
	@echo "    clean - clean all files from docs and pip routines"
	@echo "    install - install the package"
	@echo "    test - execute unit tests"
	@echo "    pictures-for-docs - plot all pictures for the documentation"
	@echo "    files-for-docs - prepare some generated files for the documentation"
	@echo "    requirements - install all requirements"

# Development environment
requirements:
	@pip install -r requirements.txt --no-cache
	@pip install -r requirements-dev.txt --no-cache
	@pip install -r docs/requirements.txt --no-cache
	@pip install -r tests/requirements.txt --no-cache

install:
	@python3 -m pip install .

clean:
	-@rm -r docs/_build
	-@rm -r docs/source/*/generated
	-@rm -r docs/source/*/*/generated
	-@rm -r docs/source/*/*/*/generated
	-@rm -r magnopy/magnopy.egg-info
	-@rm -r build
	-@rm -r dist
	-@rm -r .venv/lib/python*/site-packages/magnopy*
	-@rm -r .venv/bin/magnopy*

# Documentation and doctests
pictures-for-docs:
	@python3 dev-tools/images/origin-upstream-local.py -rd .
	@python3 dev-tools/images/positions.py -rd .

files-for-docs:
	-@rm docs/source/user-guide/cli/magnopy-lswt/help.inc
	-@rm docs/source/user-guide/cli/magnopy-optimize-sd/help.inc
	@magnopy-lswt --help > docs/source/user-guide/cli/magnopy-lswt/help.inc
	@magnopy-optimize-sd --help > docs/source/user-guide/cli/magnopy-optimize-sd/help.inc

html:
	@$(SPHINXBUILD) -M html "docs/$(SOURCEDIR)" "docs/$(BUILDDIR)" $(SPHINXOPTS)

clean-html: clean install html
	@echo "Done"

html-from-zero: clean install pictures-for-docs files-for-docs html
	@echo "Done"

# Tests
doctest:
	@$(SPHINXBUILD) -b doctest "docs/$(SOURCEDIR)" "docs/$(BUILDDIR)" $(SPHINXOPTS)

test:
	@pytest -s tests #-o log_cli=true -o log_cli_level=DEBUG
