#!/usr/bin/env bash


ls  $1 | xargs -P 6 -I NAME md5sum NAME