#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.Abstract import JavaTool

import os
from Routines.Functions import check_path


class GenotypeGVCFs(JavaTool):
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T GenotypeGVCFs", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def parse_options(self, reference, gvcf_list, output):

        options = " -nt %i" % self.threads
        options += " -R %s" % reference

        for gvcf in gvcf_list:
            options += " --variant %s" % gvcf

        options += " -o %s" % output

        return options

    def genotype(self, reference, gvcf_list, output):
        """
        java -jar GenomeAnalysisTK.jar \
           -T GenotypeGVCFs \
           -R reference.fasta \
           --variant sample1.g.vcf \
           --variant sample2.g.vcf \
           -o output.vcf
        """
        options = self.parse_options(reference, gvcf_list, output)

        self.execute(options)

