.PHONY: all html clean reset-gh-pages update-gh-pages

# builds the html documentation (run `make docs` from the repo root,
# or `make -C doc html` after `python3 makefilesegments.py`)
html:
	sphinx-build -b html -d build/doctrees . build/html

clean:
	rm -rf build

all: html

# the gh-pages targets need the orphan gh-pages branch to exist (create
# it once with a .nojekyll file, e.g.:
#   git checkout --orphan gh-pages
#   git rm -rf --cached . ; git clean -fd
#   touch .nojekyll ; git add . ; git commit -m 'init gh-pages'
#   git checkout - )
# the gh-pages worktree is build/html (avoid underscore directories,
# they can cause issues with github pages)
reset-gh-pages:
	rm -rf build/html/
	git worktree prune
	git worktree add build/html/ gh-pages
	git worktree list

update-gh-pages:
	rm -rf build/html/.buildinfo build/html/*
	$(MAKE) html
	touch build/html/.nojekyll
	git -C build/html add .
	git -C build/html commit -av -m 'update doc'
