#!/usr/bin/env bash

#echo "cat $1 | sort -k1,1 -k4,4n -k5,5"
sort -k1,1 -k4,4n -k5,5 $@
