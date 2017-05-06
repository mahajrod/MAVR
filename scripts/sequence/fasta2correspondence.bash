#!/usr/bin/env bash

#makes correspondece file with same ids

grep -P "^>" $1 | sed "s/>//;s/\s.*//;s/\(.*\)/\1\t\1/" > $2