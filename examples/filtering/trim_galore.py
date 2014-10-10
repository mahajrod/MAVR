#!/usr/bin/env python
import os
from Tools.FilterTools import trim_galore

min_length = 150
#forward_reads = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/csa/data/001/IonXpress_001.R_2013_02_25_18_36_15_user_SN1-9-Alsu_yeast_DNA.fastq"
forward_reads = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/csa/data/002/IonXpress_002.R_2013_02_25_18_36_15_user_SN1-9-Alsu_yeast_DNA.fastq"
forward_trim = 10          # 5' trim
quality_score = "phred33"  # phred33 or phred64
adapter = "ATCACCGACTGCCCATAGAGAGGAAAGCGGAGGCGTAGTGG"  # adapter sequence to be trimmed

#adapter = "GGTGATGCGGAGGCG"
quality_treshold = 20
output_folder = "./002"       # output to current folder

if output_folder != "./":
    os.mkdir(output_folder)

trim_galore(min_length, forward_reads, forward_trim,
                quality_score=quality_score,
                adapter=adapter,
                quality_treshold=quality_treshold,
                output_folder=output_folder)