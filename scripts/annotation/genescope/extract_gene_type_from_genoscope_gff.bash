#!/usr/bin/env bash

grep -v "^#" $1 | grep -P "\tgene\t" | sed 's/.*\t//' | sed 's/.*;gene_type=\([^;]\+\).*gene_name=\([^;]\+\).*/\2\t\1/'