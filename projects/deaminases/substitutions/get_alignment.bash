#!/usr/bin/env bash
WORKDIR="/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/polymorphisms/"

cd ${WORKDIR}

mkdir - p alignments/

for FILE in `ls proteins/`;
    do

    mafft --maxiterate 1000 --globalpair proteins/${FILE} > alignment/${FILE}

    done
