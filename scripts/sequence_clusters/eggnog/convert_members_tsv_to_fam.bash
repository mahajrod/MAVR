#!/usr/bin/env bash

awk -F'\t' '{printf "%s\t%s\n", $2, $6}' $1 > $2
