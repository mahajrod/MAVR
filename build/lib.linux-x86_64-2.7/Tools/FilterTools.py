#!/usr/bin/env python
import os


def trim_galore(min_length,
                forward_reads,
                forward_trim,
                reverse_reads=None,
                reverse_trim=None,
                quality_score="phred33",
                adapter="AGATCGGAAGAGC",
                quality_treshold=20,
                output_folder="trimmed"):
    # if forward_trim == None skip 5' trimming
    trim_str = ""

    if forward_trim:
        trim_str = "--clip_R1 %i" % forward_trim
        if reverse_reads:
            trim_str += " --clip_R2 %i" % reverse_trim

    reads = forward_reads
    if reverse_reads:
        reads += " %s" % reverse_reads
        reads = "--paired" + reads

    os.system("trim_galore -a %s  --length %i --%s --dont_gzip %s -q %i -o %s %s"
              % (adapter, min_length, quality_score, trim_str, quality_treshold, output_folder, reads))

    os.chdir("%s" % output_folder)
    os.system("fastqc -t 2 --nogroup *.f*q")
    os.system("cd ..")