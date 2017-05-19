#!/usr/bin/env bash

samtools view $1 | sed  's/\(^[^\t]\+\).*cx:i:\([A-Za-z0-9]\+\).*/\1\t\2/' | less