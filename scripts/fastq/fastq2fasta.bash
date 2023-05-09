#!/usr/bin/env bash

cat $1 | paste - - - - | awk '{printf ">%s\n%s\n", $1, $2}'
