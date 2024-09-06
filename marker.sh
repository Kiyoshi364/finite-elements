#!/bin/sh

set -xe

filenum=${1:-00}
file="notas/${filenum}-notas.md"

marker --preview "${file}" &
