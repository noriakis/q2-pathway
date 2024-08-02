.PHONY: all viz-kegg-pathway install dev clean distclean test

PYTHON ?= python

all: viz-kegg-pathway

q2_pathway/assets/dist:
	cd q2_pathway/assets && \
	npm install --no-save && \
	npm run build && \
        cp licenses/* dist/

viz-kegg-pathway: q2_pathway/assets/dist

install: all
	$(PYTHON) setup.py install

dev: all
	pip install -e .

clean: distclean
	rm -rf q2_pathway/assets/node_modules

distclean:
	rm -rf q2_pathway/assets/dist

test: all
	py.test