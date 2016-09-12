#!/usr/bin/env python
import os
from Tools.Abstract import JavaTool


class CombineVariants(JavaTool):

    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T CombineVariants", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html

    def combine_from_same_source(self, reference_file, input_filelist, output_vcf):

        options = " -R %s" % reference_file
        options += " -nt %i" % self.threads
        options += " --variant %s" % (input_filelist if isinstance(input_filelist, str) else " -variant ".join(input_filelist))
        options += " -o %s" % output_vcf

        self.execute(options=options)

    def combine_from_different_sources(self, gatk_dir, reference, output_file,
                                       file_list, merge_option="REQUIRE_UNIQUE"):
        variants_string = " ".join(list(map(lambda x: "--variant %s" % x, file_list)))

        os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T CombineVariants %s -o %s -genotypeMergeOptions %s"
                  % (gatk_dir, reference, variants_string, output_file, merge_option))
