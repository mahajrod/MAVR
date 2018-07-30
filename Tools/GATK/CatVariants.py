#!/usr/bin/env python

__author__ = 'mahajrod'
import shutil
import numpy as np
from Tools.Abstract import JavaTool


class CatVariants(JavaTool):
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def parse_options(self, reference, gvcf_list, output, input_is_sorted=False, extension_list=["g.vcf",]):

        #options = " -nt %i" % self.threads # bugs in tool - fails in multithreading mode
        options = " -R %s" % reference

        #gvcf_file_list = self.make_list_of_path_to_files(gvcf_list)
        #print "GGGGGGGGGGGGGGGGGGG"
        #print extension_list
        gvcf_file_list = self.make_list_of_path_to_files_by_extension(gvcf_list, extension_list=extension_list,
                                                                      recursive=False, return_absolute_paths=True)

        for gvcf in gvcf_file_list:
            options += " -V %s" % gvcf

        options += " -out %s" % output
        options += " --assumeSorted" if input_is_sorted else ""

        return options

    def combine_gvcf(self, reference, gvcf_list, output, input_is_sorted=False, extension_list=["g.vcf",],
                     tmp_dir="./tmp_combine_gvcf/", max_files_per_merging=50, iteration=0, threads=None,
                     remove_intermediate_files=False):
        """
        java -jar GenomeAnalysisTK.jar \
           -T GenotypeGVCFs \
           -R reference.fasta \
           --variant sample1.g.vcf \
           --variant sample2.g.vcf \
           -o output.vcf
        """
        #print "YYYYYYYYYYYYYYYYYYYYYY"
        #print extension_list
        filtered_gvcf_list = []
        for filename in gvcf_list:
            for extension in extension_list:
                if extension == filename[-len(extension):]:
                    filtered_gvcf_list.append(filename)
                    break
                
        if len(filtered_gvcf_list) <= max_files_per_merging:
            options = self.parse_options(reference, filtered_gvcf_list, output, input_is_sorted, extension_list=extension_list)
            self.execute(options, runtype="cp")
            if remove_intermediate_files:
                shutil.rmtree(tmp_dir, ignore_errors=True)

        else:
            self.safe_mkdir(tmp_dir)
            iteration_dir = "%s/iteration_%i/" % (tmp_dir, iteration)
            self.safe_mkdir(iteration_dir)

            number_of_files = len(filtered_gvcf_list)

            bins = np.arange(0, number_of_files, max_files_per_merging)
            #print(bins)
            if bins[-1] != number_of_files:
                if number_of_files - bins[-1] < 2:
                    bins[-1] = number_of_files
                else:
                    bins = np.append(bins, number_of_files)

            output_file_list = []
            options_list = []

            merged_files = 0
            for i in range(0, len(bins)-1):
                output_file = "%s/%i.g.vcf" % (iteration_dir, i)
                output_file_list.append(output_file)
                #print(bins[i], bins[i+1])

                merged_files += bins[i+1] - bins[i]
                options_list.append(self.parse_options(reference,
                                                       filtered_gvcf_list[bins[i]:bins[i+1]],
                                                       output_file,
                                                       input_is_sorted, extension_list=extension_list))
            print("%i/%i files will be merged" % (merged_files, number_of_files))

            self.parallel_execute(options_list, threads=threads, runtype="cp")

            self.combine_gvcf(reference, output_file_list, output, input_is_sorted=input_is_sorted,
                              extension_list=extension_list,
                              tmp_dir=tmp_dir,
                              max_files_per_merging=max_files_per_merging, iteration=iteration+1)


