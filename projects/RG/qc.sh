#!/usr/bin/env bash

#for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123; do  mkdir -p ../../unpacked/GR00/${SAMPLE}; array=($(ls ${SAMPLE}/ | sed 's/.gz//')); echo ${array[0]};  gunzip -c ${SAMPLE}/${array[0]}.gz > ../../unpacked/GR00/${SAMPLE}/${array[0]} &  echo ${array[1]}; gunzip -c ${SAMPLE}/${array[1]}.gz > ../../unpacked/GR00/${SAMPLE}/${array[1]};  done

for SAMPLE in GR0076 GR0082 GR0081 GR0122 GR0121 GR0123;
    do
    mkdir -p ../../unpacked/GR00/${SAMPLE};
    FILES=($(ls ${SAMPLE}/ | sed 's/.gz//'));
    echo ${FILES[0]};
    gunzip -c ${SAMPLE}/${FILES[0]}.gz > ../../unpacked/GR00/${SAMPLE}/${FILES[0]} &
    echo ${FILES[1]};
    gunzip -c ${SAMPLE}/${FILES[1]}.gz > ../../unpacked/GR00/${SAMPLE}/${FILES[1]};
    done

