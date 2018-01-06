# Created by G. Peter Lepage (Cornell University) on 2016-08-20.
# Copyright (c) 2016-17 G. Peter Lepage.
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
VERSION = `python -c 'import g2tools; print g2tools.__version__'`

install :
	$(PIP) install . --user

install-sys :
	$(PIP) install .

uninstall :			# mostly works (may leave some empty directories)
	$(PIP) uninstall g2tools

try:
	$(PYTHON) setup.py install --user --record files-g2tools.$(PYTHONVERSION)

untry:
	- cat files-g2tools.$(PYTHONVERSION) | xargs rm -rf

doc-html:
	rm -rf doc-build/html; sphinx-build -b html doc-build/source doc-build/html

doc-pdf:
	rm -rf doc-build/g2tools.pdf
	sphinx-build -b latex doc-build/source doc-build/latex
	cd doc-build/latex; make g2tools.pdf; mv g2tools.pdf ..

doc-zip doc.zip:
	cd doc-build/html; zip -r doc *; mv doc.zip ../..

doc-all: doc-html doc-pdf doc-zip

sdist:			# source distribution
	$(PYTHON) setup.py sdist

.PHONY: tests

tests test-all:
	$(MAKE) -C tests PYTHON=$(PYTHON) tests

run run-examples:
	$(MAKE) -C examples PYTHON=$(PYTHON) run

register-pypi:
	python setup.py register # use only once, first time

upload-pypi:
	python setup.py sdist upload

upload-git:
	echo  "version $(VERSION)"
	make doc-all
	git commit -a -m "prep documentation for upload"
	git push origin master

tag-git:
	echo  "version $(VERSION)"
	git tag -a v$(VERSION) -m "version $(VERSION)"
	git push origin v$(VERSION)

test-download:
	-$(PIP) uninstall g2tools
	$(PIP) install g2tools --no-cache-dir

clean:
	rm -f -r build
	rm -rf __pycache__
	rm -f *.so *.tmp *.pyc *.prof *.c .coverage doc.zip
	rm -f -r dist
	rm -f -r doc-build/build
	rm -f -r src/g2tools/*.c
	$(MAKE) -C doc-build/source clean
	$(MAKE) -C tests clean
	$(MAKE) -C examples clean


