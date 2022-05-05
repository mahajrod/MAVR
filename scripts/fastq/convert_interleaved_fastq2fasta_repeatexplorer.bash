#!/usr/bin/env bash

sed 's/^@/>/;s/ .*//' $1  | awk 'NR%8==1 {printf "%s/1\n", $0}; NR%8==2 {print $0}; NR%8==5 {printf "%s/2\n", $0}; NR%8==6 {print $0} '