#!/bin/sh

set -xe

DEFAULT_FILE="main"
file="${1:-$DEFAULT_FILE}"

mkdir -p folder-build
pushd folder-build
typst compile ../"$file".typ "$file".pdf
if [[ "$file" = "$DEFAULT_FILE" ]] then
  pdftotext "$file".pdf "$file".txt
  mv "$DEFAULT_FILE".pdf ../o-report.pdf
  cat "$file".txt \
    | aspell -d pt_BR \
    --add-wordlists=../wordlist.txt \
    --ignore 2 \
    list \

fi

popd
