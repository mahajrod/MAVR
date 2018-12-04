#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
from Tools.Samtools import SamtoolsV1
from Tools.Abstract import JavaTool


class MarkDuplicates(JavaTool):

    def __init__(self, java_path="", max_threads=4, jar_path="", max_memory="1g"):
        jar = "picard.jar MarkDuplicates"
        JavaTool.__init__(self, jar, java_path=java_path, max_threads=max_threads,
                          jar_path=jar_path, max_memory=max_memory)

    def parse_options(self, input_bam, output_bam, stat_file):

        options = " I= %s" % input_bam
        options += " O= %s" % output_bam
        options += " M= %s" % stat_file

        return options

    def run(self, input_bam, output_bam, stat_file):

        options = self.parse_options(input_bam, output_bam, stat_file)

        self.execute(options)

    def mkdup(self, input_bam, output_prefix):
        output_bam = "%s.bam" % output_prefix
        stat_file = "%s.stat" % output_prefix

        options = self.parse_options(input_bam, output_bam, stat_file)

        self.execute(options=options)
        SamtoolsV1.index(output_bam)

if __name__ == "__main__":
    pass