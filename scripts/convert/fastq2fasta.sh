#!/usr/bin/env bash

#INPUT_FILES=$1
#replace comma by space
#INPUT_FILES=${INPUT_FILES//,/ }

#cat ${INPUT_FILES} | awk '{if(NR%4==1 || NR%4==2) print}' | sed 's/^@/\>/'

cat $@ | awk '{if(NR%4==1 || NR%4==2) print $0 }' | sed 's/^@/\>/'