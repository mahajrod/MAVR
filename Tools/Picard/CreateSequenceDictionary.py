#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import JavaTool


class CreateSequenceDictionary(JavaTool):

    def __init__(self, java_path="", max_threads=4, jar_path="", max_memory="1g"):
        jar = "CreateSequenceDictionary.jar"
        JavaTool.__init__(self, jar, java_path=java_path, max_threads=max_threads,
                          jar_path=jar_path, max_memory=max_memory)

    def make_fasta_dict(self, fasta_file, dict_name):

        options = " R= %s" % fasta_file
        options += " O= %s" % dict_name

        self.execute(options)

if __name__ == "__main__":
    pass
