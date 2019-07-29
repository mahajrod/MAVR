#!/usr/bin/env bash

cat $0 | sort -k1,1 -k4,4n -k5,5
