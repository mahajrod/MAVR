#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.Abstract import JavaTool

import os
from Routines.Functions import check_path


class CatVariants(JavaTool):
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def parse_options(self, reference, gvcf_list, output, input_is_sorted=False, extension_list=["g.vcf",]):

        #options = " -nt %i" % self.threads # bugs in tool - fails in multithreading mode
        options = " -R %s" % reference

        #gvcf_file_list = self.make_list_of_path_to_files(gvcf_list)

        gvcf_file_list = self.make_list_of_path_to_files_by_extension(gvcf_list, extension_list=extension_list,
                                                                      recursive=False, return_absolute_paths=True)

        for gvcf in gvcf_file_list:
            options += " -V %s" % gvcf

        options += " -out %s" % output
        options += " --assumeSorted" if input_is_sorted else ""

        return options

    def combine_gvcf(self, reference, gvcf_list, output, input_is_sorted=False, extension_list=["g.vcf",]):
        """
        java -jar GenomeAnalysisTK.jar \
           -T GenotypeGVCFs \
           -R reference.fasta \
           --variant sample1.g.vcf \
           --variant sample2.g.vcf \
           -o output.vcf
        """
        options = self.parse_options(reference, gvcf_list, output, input_is_sorted, extension_list=extension_list)

        self.execute(options, runtype="cp")

