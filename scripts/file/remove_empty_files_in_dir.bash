#!/usr/bin/env bash


for FILE in $1/*; do FSIZE=`ls -l ${FILE} | awk '{print $5}'`;if [ ${FSIZE} -eq 0 ]; then echo ${FILE}; rm ${FILE}; fi; done
