#!/usr/bin/env python
import os
from Tools.Abstract import Tool


class Bowtie2(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bowtie2", path=path, max_threads=max_threads)

    def index(self, reference, index_name):
        options = "%s %s" % (reference, index_name)
        self.execute(options, cmd="bowtie2-build")

    def align(self,
              bowtie2_index,
              right_reads_list=None,
              left_reads_list=None,
              unpaired_reads_list=None,
              quality_score="phred33",
              alignment_mode="very-sensitive",
              output_file="alignment.sam",
              find_discordant_alignments=True,
              find_separated_alignments=True):

        options = " -p %i" % self.threads
        options += " --%s" % alignment_mode
        options += " --%s" % quality_score
        options += " -x %s" % bowtie2_index
        options += " --no-discordant" if not find_discordant_alignments else ""
        options += " --no-mixed" if not find_separated_alignments else ""
        options += " -1 %s -2 %s" % (",".join(right_reads_list), ",".join(left_reads_list)) \
            if right_reads_list and left_reads_list else ""
        options += " -U %s" % ",".join(unpaired_reads_list) if unpaired_reads_list else ""
        options += " > %s" % output_file

        self.execute(options)

