#!/usr/bin/env bash

grep -P '\tCDS\t' $1 | grep 'protein_id=' | sed 's/.*gene=\([A-Za-z0-9\.\-\_]\+\).*protein_id=\([A-Za-z0-9\.\-\_]\+\).*/\1\t\2/'  | sort | uniq