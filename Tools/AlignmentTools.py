#!/usr/bin/env python
import os


class Bowtie2():

    def align(self,
              bowtie2_index,
              forward_reads,
              reverse_reads=None,
              max_threads=5,
              quality_score="phred33",
              alignment_mode="very-sensitive",
              output_file="alignment.sam"):

        if reverse_reads:
            reads = "-1 %s -2 %s" % (forward_reads, reverse_reads)
        else:
            reads = "-U %s" % forward_reads

        os.system("bowtie2 --%s --%s -p %i -x %s %s > %s"
                      % (alignment_mode, quality_score, max_threads, bowtie2_index, reads, output_file))