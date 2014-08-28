#!/usr/bin/env python
import os
from Tools.Abstract import Tool


class Bowtie2(Tool):
    def __init__(self, path="", max_threads=4):
        super(Tool).__init__("bowtie2", path=path, max_threads=max_threads)

    def index(self, reference, index_name):
        options = "%s %s" % (reference, index_name)
        self.execute(options, cmd="bowtie2-build")
        #os.system("bowtie2-build %s %s" % (reference, index_name))

    def align(self,
              bowtie2_index,
              forward_reads,
              reverse_reads=None,
              quality_score="phred33",
              alignment_mode="very-sensitive",
              output_file="alignment.sam",
              find_discordant_alignments=True,
              find_separated_alignments=True):

        if reverse_reads:
            reads = "-1 %s -2 %s" % (forward_reads, reverse_reads)
        else:
            reads = "-U %s" % forward_reads

        options = ""
        if not find_discordant_alignments:
            options += " --no-discordant"
        if not find_separated_alignments:
            options += " --no-mixed"

        options = "--%s --%s -p %i %s -x %s %s > %s" %\
                  (alignment_mode, quality_score, self.threads, options, bowtie2_index, reads, output_file)

        self.execute(options)

