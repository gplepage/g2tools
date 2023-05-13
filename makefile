# Created by G. Peter Lepage (Cornell University) on 2016-08-20.
# Copyright (c) 2016-23 G. Peter Lepage.
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
SRCFILES := $(shell ls setup.py src/g2tools/*.py)

install-user :
	$(PIP) install . --user --no-cache-dir

install install-sys :
	$(PIP) install . --no-cache-dir

uninstall :			# mostly works (may leave some empty directories)
	$(PIP) uninstall g2tools

update: 
	make uninstall install

.PHONY : doc

doc-html doc: 
	make doc/html/index.html

doc/html/index.html : $(DOCFILES) $(SRCFILES) setup.cfg
	sphinx-build -b html doc/source doc/html

clear-doc:
	rm -rf doc/html; 

sdist:			# source distribution
	$(PYTHON) -m build --sdist

.PHONY: tests

tests test-all:
	$(MAKE) -C tests PYTHON=$(PYTHON) tests

run run-examples:
	$(MAKE) -C examples PYTHON=$(PYTHON) run

upload-twine:
	twine upload dist/g2tools-$(VERSION).tar.gz

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


