#!/usr/bin/env python
import os
from Tools.Abstract import JavaTool


class ValidateVariants(JavaTool):
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_filters_VariantFiltration.html
    # default filters for indel and snp filtration were taken from GATK BestPractice
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T ValidateVariants", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    @staticmethod
    def parse_common_options(reference, input_vcf, exclude_list=[], dbsnp=None, input_is_gvcf=False):
        options = " -R %s" % reference
        options += " -V %s" % input_vcf
        options += " --dbsnp %s" % dbsnp if dbsnp else ""
        for entry in exclude_list:
            options += " ----validationTypeToExclude %s" % entry

        options += " --validateGVCF" if input_is_gvcf else ""

        return options

    def test_vcf_format(self, reference, input_vcf, input_is_gvcf=False):

        options = self.parse_common_options(reference, input_vcf, exclude_list=["ALL", ], input_is_gvcf=input_is_gvcf)

        self.execute(options=options)

    def index_vcf(self, reference, input_vcf, input_is_gvcf=False):
        self.test_vcf_format(reference, input_vcf, input_is_gvcf=input_is_gvcf)

