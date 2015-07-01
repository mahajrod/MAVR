#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import JavaTool


class AddOrReplaceReadGroups(JavaTool):

    def __init__(self, java_path="", max_threads=4, jar_path="", max_memory="1g"):
        jar = "AddOrReplaceReadGroups.jar"
        JavaTool.__init__(self, jar, java_path=java_path, max_threads=max_threads,
                          jar_path=jar_path, max_memory=max_memory)

    def add_read_groups(self, input_bam, output_bam, RGID, RGLB, RGPL, RGSM, RGPU):

        options = " I= %s" % input_bam
        options += " O= %s" % output_bam
        options += " SORT_ORDER=coordinate"
        options += " RGID=%s" % RGID
        options += " RGLB=%s" % RGLB
        options += " RGPL=%s" % RGPL
        options += " RGSM=%s" % RGSM
        options += " RGPU=%s" % RGPU
        options += " CREATE_INDEX=True"

        self.execute(options)

if __name__ == "__main__":
    pass