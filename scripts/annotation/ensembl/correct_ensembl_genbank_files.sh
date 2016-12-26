#!/usr/bin/env bash

#Tested on ferret file from Ensembl 87

#ACCESSION should be same to VERSION and LOCUS id due to requirenments of biopython genbank parser
#this script edits:
#   VERSION field - added point before last digit
#   LOCUS

#Correction from


# LOCUS       GL896898.1 52375790 bp DNA HTG 24-NOV-2016
# DEFINITION  Mustela putorius furo scaffold GL896898.1 MusPutFur1.0 full sequence 1..52375790
#             reannotated via EnsEMBL
# ACCESSION   scaffold:MusPutFur1.0:GL896898.1:1:52375790:1
# VERSION     GL8968981
#
# to
#
# LOCUS       GL896898.1 52375790 bp DNA HTG 24-NOV-2016
# DEFINITION  Mustela putorius furo scaffold GL896898.1 MusPutFur1.0 full sequence 1..52375790
#             reannotated via EnsEMBL
# ACCESSION   GL896898.1
# VERSION     GL896898.1


sed 's/^\(ACCESSION \+\)scaffold:[a-zA-Z0-9.]\+:\([a-zA-Z0-9.]\+\):[0-9]\+:[0-9]\+:[0-9]\+/\1\2/;s/^\(VERSION \+[a-zA-Z0-9]\+\)\([0-9]\)/\1.\2/' $1
