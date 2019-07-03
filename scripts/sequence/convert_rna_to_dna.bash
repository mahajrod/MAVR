#!/usr/bin/env bash

sed "/^[^>]/s/U/T/g;/^[^>]/s/u/t/g" $1
