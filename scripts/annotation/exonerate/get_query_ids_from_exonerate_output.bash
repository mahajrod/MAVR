#!/usr/bin/env bash

grep -h "         Query:" $1 | sed -r 's/\s+Query:\s+//; s/\s+.*//' | sort | uniq
