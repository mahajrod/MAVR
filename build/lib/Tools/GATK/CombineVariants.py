#!/usr/bin/env python
import os


class CombineVariants():
    # TODO: rewrite more clear
    # http://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_sting_gatk_walkers_variantutils_CombineVariants.html

    def combine_from_same_source(self, gatk_dir, reference, input_filelist, output_file, num_threads=4):
        in_files = ""
        for entry in input_filelist:
            in_files += " --variant %s" % entry
        os.system("java -jar %sGenomeAnalysisTK.jar -R %s -T CombineVariants -nt %i %s -o %s"
                  % (gatk_dir, reference, num_threads, in_files, output_file))

    def combine_from_different_sources(self, gatk_dir, reference, output_file,
                                       file_list, merge_option="REQUIRE_UNIQUE"):
        variants_string = " ".join(list(map(lambda x: "--variant %s" % x, file_list)))

        os.system("java -Xmx2g -jar %sGenomeAnalysisTK.jar -R %s -T CombineVariants %s -o %s -genotypeMergeOptions %s"
                  % (gatk_dir, reference, variants_string, output_file, merge_option))
