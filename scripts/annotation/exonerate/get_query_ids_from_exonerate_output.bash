#!/usr/bin/env bash

grep -h "         Query:" splited_output/* | sed -r 's/\s+Query:\s+//; s/\s+.*//' | sort | uniq
