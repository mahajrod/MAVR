#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import shutil

from Tools.Abstract import Tool

from Parsers.TRF import CollectionTRF
from Routines.File import split_filename, save_mkdir
from Tools.LinuxTools import CGAS


class FastqDump(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "fastq-dump", path=path, max_threads=max_threads)

    def unpack(self, input_sra, output_dir, paired_input=True, retain_original_ids=True):

        options = " --outdir %s" % output_dir
        options += " --split-3" if paired_input else ""
        options += " --origfmt" if retain_original_ids else ""
        options += " %s" % input_sra

        self.execute(options=options)

    def parallel_unpack(self, input_dir, output_dir, sra_id_list=None, paired_input=True, retain_original_ids=True):

        sra_id_list = sra_id_list if sra_id_list else os.listdir()

        common_options = " --split-3" if paired_input else ""
        common_options += " --origfmt" if retain_original_ids else ""

        options_list = []

        for sra_id in sra_id_list:
            out_dir = "%s/%s/" % (output_dir, sra_id)
            input_sra = "%s/%s/%s.sra" % (input_dir, sra_id, sra_id)

            self.safe_mkdir(out_dir)

            options = common_options
            options += " --outdir %s" % out_dir
            options += " %s" % input_sra

            options_list.append(options)

        self.parallel_execute(options_list=options_list)

if __name__ == "__main__":
    pass
