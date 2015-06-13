#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Spades(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "spades.py", path=path, max_threads=max_threads)

    def assembly(self, pacbio_reads=None, sanger_reads=None, only_error_correction=False,
                 only_assembler=False, careful=True, single_cell=False, iontorrent=False,
                 trusted_contigs=None, untrusted_contigs=None, out_dir="spades"):
        options = " -t %i" % self.threads
        options += " --sanger %s" % (sanger_reads) if sanger_reads else ""
        options += " --pacbio %s" % (pacbio_reads) if pacbio_reads else ""
        options += " --only-error-correction" if only_error_correction else ""
        options += " --only_assembler" if only_assembler else ""
        options += " --careful" if careful else ""
        options += " --sc" if single_cell else ""
        options += " --iontorrent" if iontorrent else ""
        options += " -o %s" % out_dir
        options += " --trusted-contigs %s" % trusted_contigs if trusted_contigs else ""
        options += " --untrusted-contigs %s" % untrusted_contigs if untrusted_contigs else ""
        options += " "
        options += " "
        options += " "
        self.execute(options)

    def error_correction(self):
        pass

if __name__ == "__main__":
    pass