PYTHON = python

run run-all run-examples:
	-$(PYTHON) amu-example.py > amu-example.tmp
	-diff -w -I 't=' -I svdcut amu-example.tmp amu-example.out

update-all:
	-$(PYTHON) amu-example.py > amu-example.out

clean:
	rm -f *.pyc *.p *.tmp *.prof