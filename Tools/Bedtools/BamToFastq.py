#!/usr/bin/env python
from Tools.Abstract import Tool


class BamToFastq(Tool):
    """
    T
    """

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bedtools bamtofastq", path=path, max_threads=max_threads)

    def convert(self, in_bam, out_left_fastq, out_right_fastq=None, use_bam_mate_tags=False):

        options = " -tags" if use_bam_mate_tags else ""
        options += " -i %s" % in_bam
        options += " -fq %s" % out_left_fastq
        options += " -fq2 %s" % out_right_fastq if out_right_fastq else ""

        self.execute(options)

