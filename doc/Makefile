.PHONY: all html clean

# builds the html documentation (run `make docs` from the repo root,
# or `make -C doc html` after `python3 makefilesegments.py`)
html:
	sphinx-build -b html -d _build/doctrees . _build/html

clean:
	rm -rf _build

all: html
