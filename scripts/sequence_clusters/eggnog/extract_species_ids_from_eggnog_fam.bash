#!/usr/bin/env bash

awk -F '\t' '{print $2}' $1 | tr ',' '\n' | sed 's/\..*//' | sort | uniq
