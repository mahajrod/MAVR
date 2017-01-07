#!/usr/bin/env bash

samtools view $1 | awk  '{print $1}' |  awk 'BEGIN {B=1}; { if (B == 2) {if($1 == FORWARD_NAME) {B=1;} else {print $1} } else {FORWARD_NAME=$1; B=2};   }'