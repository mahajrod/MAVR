#!/usr/bin/env bash

grep gene_biotype=protein_coding $1 | grep -P "\tgene\t" | grep -vP "\tGnomon\t" |cut -f 9 | cut -d ";" -f 1 | cut -d = -f 2