#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
from Tools.Abstract import JavaTool


class CreateSequenceDictionary(JavaTool):

    def __init__(self, java_path="", max_threads=4, jar_path="", max_memory="1g"):
        #JavaTool.__init__(self, "CreateSequenceDictionary.jar", java_path=java_path, max_threads=max_threads,
        #                  jar_path=jar_path, max_memory=max_memory)
        JavaTool.__init__(self, "picard.jar CreateSequenceDictionary", java_path=java_path, max_threads=max_threads,
                          jar_path=jar_path, max_memory=max_memory)

    def make_fasta_dict(self, fasta_file, dict_name):

        options = " R= %s" % fasta_file
        options += " O= %s" % dict_name

        self.execute(options)

    def check_for_fasta_dict(self, fasta_file):

        fasta_path_list = self.split_filename(fasta_file)
        dict_file = "%s/%s.dict" % (fasta_path_list[0], fasta_path_list[1])

        if not os.path.exists(dict_file):
            self.make_fasta_dict(fasta_file, dict_file)


if __name__ == "__main__":
    pass
