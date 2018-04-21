#!/usr/bin/env bash

sed 's/\ /_/g;s/\(.*\)\t\(.*\)/s\/>\1\.\/>\2.\//' $1

#sed 's/\(.*\)\t\(.*\)/s\/>\1\.\/>\2.\//' species.syn  > species.sed
