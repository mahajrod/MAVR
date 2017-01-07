#!/usr/bin/env bash


sed 's/^\(>.*\)@.*/\1/' $1 > $2
