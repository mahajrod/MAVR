#!/usr/bin/env bash

sed 's/\(.*\)\t\(.*\)/s\/\1\.\/\2.\//' $1

#sed 's/\(.*\)\t\(.*\)/s\/>\1\.\/>\2.\//' species.syn  > species.sed
