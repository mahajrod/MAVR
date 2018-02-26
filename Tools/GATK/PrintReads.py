#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.GATK.Abstract import GATKTool


class PrintReads(GATKTool):
    def __init__(self, java_path="", max_threads=4, jar_path="", max_memory="1g"):
        GATKTool.__init__(self, java_path=java_path, max_threads=max_threads, jar_path=jar_path,
                          max_memory=max_memory)

    def get_recalled_reads(self, reference, input_bam, recal_table, output_bam):
        options = " -nct %i" % self.threads
        options += " -R %s" % reference
        options += " -I %s" % input_bam
        options += " -BQSR %s" % recal_table
        options += " -o %s" % output_bam

        self.execute(options, cmd="-T PrintReads")
