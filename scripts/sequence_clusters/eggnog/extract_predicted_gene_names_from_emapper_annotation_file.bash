#!/usr/bin/env bash

awk -F'\t' '{printf "%s\t%s\n", $1, $5}' $1
