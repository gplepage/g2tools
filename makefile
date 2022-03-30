# Created by G. Peter Lepage (Cornell University) on 2016-08-20.
# Copyright (c) 2016-18 G. Peter Lepage.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# any later version (see <http://www.gnu.org/licenses/>).
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

PIP = python -m pip
PYTHON = python
PYTHONVERSION = python`python -c 'import platform; print(platform.python_version())'`
VERSION = `cd ..; python -c 'import g2tools; print(g2tools.__version__)'`

DOCFILES :=  $(shell ls doc/source/*.{rst,py})
SRCFILES := $(shell ls setup.py src/*.py)

install-user :
	$(PIP) install . --user

install install-sys :
	$(PIP) install .

uninstall :			# mostly works (may leave some empty directories)
	$(PIP) uninstall g2tools

try:
	$(PYTHON) setup.py install --user --record files-g2tools.$(PYTHONVERSION)

untry:
	- cat files-g2tools.$(PYTHONVERSION) | xargs rm -rf

doc-html:
	make doc/html/index.html

doc/html/index.html : $(DOCFILES) $(SRCFILES)
	rm -rf doc/html; sphinx-build -b html doc/source doc/html

# doc-pdf:
# 	make doc/g2tools.pdf

# doc/g2tools.pdf : $(DOCFILES) $(SRCFILES)
# 	rm -rf doc/g2tools.pdf
# 	sphinx-build -b latex doc/source doc/latex
# 	cd doc/latex; make g2tools.pdf; mv g2tools.pdf ..

doc-zip doc.zip:
	cd doc/html; zip -r doc *; mv doc.zip ../..

doc-all: doc-html # doc-pdf

sdist:			# source distribution
	$(PYTHON) setup.py sdist

.PHONY: tests

tests test-all:
	$(MAKE) -C tests PYTHON=$(PYTHON) tests

run run-examples:
	$(MAKE) -C examples PYTHON=$(PYTHON) run

register-pypi:
	python setup.py register # use only once, first time

upload-twine:
	twine upload dist/g2tools-$(VERSION).tar.gz

upload-pypi:
	python setup.py sdist upload

upload-git:
	echo  "version $(VERSION)"
	make doc-html # doc-pdf
	git diff --exit-code
	git diff --cached --exit-code
	git push origin master

tag-git:
	echo  "version $(VERSION)"
	git tag -a v$(VERSION) -m "version $(VERSION)"
	git push origin v$(VERSION)

test-download:
	-$(PIP) uninstall g2tools
	$(PIP) install g2tools --no-cache-dir

test-readme:
	python setup.py --long-description | rst2html.py > README.html

clean:
	rm -f -r build
	rm -rf __pycache__
	rm -f *.so *.tmp *.pyc *.prof *.c .coverage doc.zip
	rm -f -r dist
	rm -f -r doc/build
	rm -f -r src/g2tools/*.c
	$(MAKE) -C doc/source clean
	$(MAKE) -C tests clean
	$(MAKE) -C examples clean


