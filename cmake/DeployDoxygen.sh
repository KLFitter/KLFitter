#!/bin/sh
#
# Required variables:
# - TRAVIS_BUILD_NUMBER : The number of the current build.
# - TRAVIS_COMMIT       : The commit that the current build is testing.
# - DOC_DOXYFILE        : The Doxygen configuration file.
# - DOC_REPO            : The documentation repository.
# - DOC_TOKEN           : Secure token to the github repository.
#

set -e
mkdir doc_tmp && cd doc_tmp

git clone https://git@github.com/KLFitter/$DOC_REPO
cd $DOC_REPO

git config user.name "Travis CI"
git config user.email "travis@travis-ci.org"

rm -rf *
echo "" > .nojekyll
doxygen $DOXYFILE 2>&1 | tee doxygen.log

if [ -d "html" ] && [ -f "html/index.html" ]; then

    echo 'Uploading documentation to the repository ...'
    git add --all
    git commit -m "Deploy documentation; travis build: ${TRAVIS_BUILD_NUMBER}" -m "Commit: ${TRAVIS_COMMIT}"
    git push --force "https://${DOC_TOKEN}@github.com/KLFitter/${DOC_REPO}" > /dev/null 2>&1
else
    echo '' >&2
    echo 'Warning: No documentation (html) files have been found!' >&2
    echo 'Warning: Not going to push the documentation to GitHub!' >&2
    exit 1
fi
