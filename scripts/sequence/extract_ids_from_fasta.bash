#!/usr/bin/env bash

grep -P "^>" $1 | sed 's/^>//;s/ .*//'
