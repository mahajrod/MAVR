#!/usr/bin/env python
import os

from Tools.Abstract import Tool
from CustomCollections.GeneralCollections import IdSet


class RepeatModeler(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "RepeatModeler", path=path, max_threads=max_threads)

    def build_db(self, database_name, fasta_dir=None, file_with_filenames=None, engine="ncbi"):

        if (not fasta_dir) and (not file_with_filenames):
            raise ValueError("Neither fasta file/directory with fasta files or file with paths to fasta files were set")

        options = " -name %s" % database_name
        options += " -dir %s" % fasta_dir if fasta_dir else ""
        options += " -engine %s" % engine if engine else ""
        options += " -batch %s " % file_with_filenames if file_with_filenames else ""

        self.cmd(cmd="BuildDatabase", options)

    def annotate_repeats(self, database, engine="ncbi", recover_dir=None):

        options = " -pa %i" % self.threads
        options += " -database %s" % database
        options += " -engine %s" % engine if engine else ""
        options += " -recoverDir %s " % recover_dir if recover_dir else ""

        self.cmd(options)