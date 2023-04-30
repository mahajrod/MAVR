#!/usr/bin/env bash

paste <(cat $1| paste - - - -) <(cat $2 | paste - - - -)  | tr "\t" "\n"
