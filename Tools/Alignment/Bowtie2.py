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
              forward_reads_list=None,
              reverse_reads_list=None,
              unpaired_reads_list=None,
              quality_score="phred33",
              alignment_mode="very-sensitive",
              find_discordant_alignments=True,
              find_separated_alignments=True,
              concordant_upper_threshold=None,
              output_prefix="alignment",
              output_format="bam",
              read_group_name="reads",
              PU="x",
              SM="sample",
              platform="Illumina",
              LB="x",
              sort_by_coordinate=True,
              sort_by_name=False,
              max_per_sorting_thread_memory=None):

        options = " -p %i" % self.threads
        options += " --%s" % alignment_mode
        options += " --%s" % quality_score
        options += " -x %s" % bowtie2_index
        options += " -X %i" % concordant_upper_threshold if concordant_upper_threshold else ""
        options += " --no-discordant" if not find_discordant_alignments else ""
        options += " --no-mixed" if not find_separated_alignments else ""
        options += " -1 %s -2 %s" % (",".join(forward_reads_list), ",".join(reverse_reads_list)) \
            if forward_right_reads_list and reverse_left_reads_list else ""
        options += " -U %s" % ",".join(unpaired_reads_list) if unpaired_reads_list else ""

        if sort_by_coordinate or sort_by_name:
            if sort_by_coordinate and sort_by_name:
                raise ValueError("Sorting by both coordinate and read name was requested")
            options += " samtools sort"
            if sort_by_name:
                options += " -n"
            options += " -@ %i" % self.threads
            options += " -m %s" % max_per_sorting_thread_memory if max_per_sorting_thread_memory else self.max_per_thread_memory
            options += " -O %s" % output_format.upper()

        options += " > %s.%s" % (output_prefix, output_format)

        self.execute(options)

