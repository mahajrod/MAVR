#!/usr/bin/env bash

sed 's/\ /_/g;s/\(.*\)\t\(.*\)/s\/\1\.\/\2.\/g/' $1

#sed 's/\(.*\)\t\(.*\)/s\/\1\.\/\2.\/g/' species.syn  > species.sed
