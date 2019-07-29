#!/usr/bin/env bash

echo "cat $0 | sort -k1,1 -k4,4n -k5,5"
cat $0 | sort -k1,1 -k4,4n -k5,5
