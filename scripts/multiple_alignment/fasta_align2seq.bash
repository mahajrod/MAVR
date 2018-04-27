#!/usr/bin/env bash

sed 's/-//g;/^\s*$/d' $1 > $2
