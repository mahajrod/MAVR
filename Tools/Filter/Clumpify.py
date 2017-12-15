#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Clumpify(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "clumpify.sh", path=path, max_threads=max_threads)

    def remove_optical_duplicates(self, forward_reads, output_prefix, reverse_reads=None, gzip_output=False,
                                  memory_limit="2g"):

        if reverse_reads:
            forward_output = "%s_1.fastq" % output_prefix + (".gz" if gzip_output else "")
            reverse_output = "%s_2.fastq" % output_prefix + (".gz" if gzip_output else "")
        else:
            output = "%s.fastq" % output_prefix + (".gz" if gzip_output else "")

        options = " in1=%s in2=%s" % (forward_reads, reverse_reads) if reverse_reads else " in=%s" % forward_reads
        options += " out1=%s out2=%s" % (forward_output, reverse_output) if reverse_reads else " out=%s" % output

        options += " dedupe"
        options += " optical"
        options += " -Xmx %s" % str(memory_limit)
        self.execute(options=options)


if __name__ == "__main__":
    pass


