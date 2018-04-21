#!/usr/bin/env bash

sed 's/\(.*\)\t\(.*\)/g;s\/\1\.\/\2.\//g' $1

#sed 's/\(.*\)\t\(.*\)/g;s\/\1\.\/\2.\//g' species.syn  > species.sed
