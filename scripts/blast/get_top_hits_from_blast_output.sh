#!/usr/bin/env bash

awk -F'\t' '{if ($1 != PREV) {print $0; PREV=$1} }' $1
