#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.Abstract import JavaTool

from Routines import VCFRoutines


class GenotypeGVCFs(JavaTool):
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog="tool_time.log"):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T GenotypeGVCFs", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def parse_options(self, reference, gvcf_list, output, extension_list=["g.vcf",]):
        options = self.parse_options_for_parallel_run(reference,
                                                      gvcf_list,
                                                      extension_list=extension_list)

        options += " -o %s" % output

        return options

    def parse_options_for_parallel_run(self, reference, gvcf_list, extension_list=["g.vcf",]):

        options = " -R %s" % reference

        for gvcf in self.make_list_of_path_to_files_by_extension(gvcf_list,
                                                                 extension_list=extension_list,
                                                                 recursive=False, return_absolute_paths=True):
            options += " --variant %s" % gvcf

        return options

    def genotype(self, reference, gvcf_list, output, extension_list=["g.vcf",]):
        """
        java -jar GenomeAnalysisTK.jar \
           -T GenotypeGVCFs \
           -R reference.fasta \
           --variant sample1.g.vcf \
           --variant sample2.g.vcf \
           -o output.vcf
        """
        options = self.parse_options(reference, gvcf_list, output, extension_list=extension_list)

        self.execute(options)

    def parallel_genotype(self, reference, gvcf_list, output_dir, output_prefix,
                          max_total_scaffold_length_per_chunk=100000,
                          max_scaffold_number_per_chunk=5, length_dict=None,
                          parsing_mode="parse", region_list=None,
                          extension_list=["g.vcf",], tmp_subdir="merging_tmp/",
                          remove_intermediate_files=False):

        self.safe_mkdir(output_dir)

        regions_list = self.prepare_region_list_by_length(max_length=max_total_scaffold_length_per_chunk,
                                                          max_seq_number=max_scaffold_number_per_chunk,
                                                          length_dict=length_dict,
                                                          reference=None if length_dict is not None else reference,
                                                          parsing_mode=parsing_mode,
                                                          output_dir="%s/regions/" % output_dir,
                                                          split_scaffolds=False) if region_list is None else region_list

        options = self.parse_options_for_parallel_run(reference, gvcf_list, extension_list=extension_list)

        output_index = 1
        options_list = []

        region_vcf_list = []

        for regions in regions_list:
            region_options = " -o %s/%s_%i.vcf" % (output_dir, output_prefix, output_index)
            region_vcf_list.append("%s/%s_%i.vcf" % (output_dir, output_prefix, output_index))
            for region in regions:
                if isinstance(region, str):
                    region_options += " -L %s" % region
                elif len(region) == 1:
                    region_options += " -L %s" % region[0]
                elif len(region) == 3:
                    region_options += " -L %s:%i-%i" % (region[0], region[1], region[2])

            options_list.append(options + region_options)
            output_index += 1

        self.parallel_execute(options_list)

        VCFRoutines.combine_same_samples_vcfs(region_vcf_list,
                                              "%s.vcf" % output_prefix,
                                              close_fd_after=False,
                                              extension_list=[".vcf", ])

        """
        CatVariants.threads = max(8, self.threads)
        CatVariants.max_memory = self.max_memory

        CatVariants.combine_gvcf(reference, output_dir, "%s.vcf" % output_prefix,
                                 input_is_sorted=True,
                                 extension_list=["vcf", ],
                                 tmp_dir="%s/%s/" %(output_dir, tmp_subdir),
                                 max_files_per_merging=100,
                                 iteration=0,
                                 threads=None,
                                 remove_intermediate_files=remove_intermediate_files)

        """