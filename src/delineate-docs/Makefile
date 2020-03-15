# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build
ARGS=$(@)

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: publish help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@echo === Executing Sphinx build ...
	@echo
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@if [ "$@" = "html" ]; then \
		echo ;\
		echo === Copying files to deploy directory ... ;\
		echo ;\
		mkdir -p "../../docs"; \
		rsync -ap --progress --delete $(wildcard $(BUILDDIR)/html/*) "../../docs/"; \
	else \
		echo ;\
		echo "(Not copying not HTML target)"; \
	fi

publish:
	git add source "../../docs"
	git commit
	git push
