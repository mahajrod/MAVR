#!/usr/bin/env python

__author__ = 'mahajrod'

from Tools.Abstract import JavaTool

from Routines import VCFRoutines

from Tools.GATK import ValidateVariants


class GenotypeGVCFs(JavaTool):
    def __init__(self,  java_path="", max_threads=4, jar_path="", max_memory=None, timelog=None):
        JavaTool.__init__(self, "GenomeAnalysisTK.jar -T GenotypeGVCFs", java_path=java_path,
                          max_threads=max_threads, jar_path=jar_path, max_memory=max_memory,
                          timelog=timelog)

    def parse_options_for_parallel_run(self, reference, gvcf_list, extension_list=["g.vcf",],
                                       disable_auto_index_creation_and_locking_when_reading_rods=True,
                                       max_alternate_alleles=None, ):

        options = " -R %s" % reference
        options += " --max_alternate_alleles %i" % max_alternate_alleles if max_alternate_alleles else ""
        options += " --disable_auto_index_creation_and_locking_when_reading_rods" if disable_auto_index_creation_and_locking_when_reading_rods else ""

        for gvcf in self.make_list_of_path_to_files_by_extension(gvcf_list,
                                                                 extension_list=extension_list,
                                                                 recursive=False, return_absolute_paths=True):
            options += " --variant %s" % gvcf
        return options

    def parse_options(self, reference, gvcf_list, output, extension_list=["g.vcf",],
                      disable_auto_index_creation_and_locking_when_reading_rods=True,
                      max_alternate_alleles=None):
        options = self.parse_options_for_parallel_run(reference,
                                                      gvcf_list,
                                                      extension_list=extension_list,
                                                      disable_auto_index_creation_and_locking_when_reading_rods=disable_auto_index_creation_and_locking_when_reading_rods,
                                                      max_alternate_alleles=max_alternate_alleles)

        options += " -o %s" % output

        return options

    def genotype(self,
                 reference,
                 gvcf_list,
                 output_prefix,
                 extension_list=["g.vcf",],
                 handling_mode="local",
                 max_memory_per_node=None,
                 job_name=None,
                 log_prefix=None,
                 error_log_prefix=None,
                 modules_list=None,
                 environment_variables_dict=None,
                 max_alternate_alleles=None,
                 max_running_time=None):
        """
        java -jar GenomeAnalysisTK.jar \
           -T GenotypeGVCFs \
           -R reference.fasta \
           --variant sample1.g.vcf \
           --variant sample2.g.vcf \
           -o output.vcf
        """
        output = "%s.vcf" % output_prefix
        options = self.parse_options(reference, gvcf_list, output, extension_list=extension_list,
                                     max_alternate_alleles=max_alternate_alleles)

        if handling_mode == 'local':
            self.execute(options)
        elif handling_mode == "slurm":
            self.timelog = None
            slurm_cmd = self.execute(options, capture_output=False, runtype="jar", generate_cmd_string_only=True)

            last_job_id = self.slurm_run_job(job_name,
                                             log_prefix,
                                             slurm_cmd,
                                             error_log_prefix,
                                             "%s.slurm" % output_prefix,
                                             modules_list=modules_list,
                                             max_running_time=max_running_time,
                                             environment_variables_dict=environment_variables_dict,
                                             max_memory_per_node=max_memory_per_node)

            return last_job_id

    def parallel_genotype(self, reference, gvcf_list, splited_dir, splited_prefix, output_vcf,
                          max_total_scaffold_length_per_chunk=100000,
                          max_scaffold_number_per_chunk=5, length_dict=None,
                          parsing_mode="parse", region_list=None,
                          extension_list=["g.vcf",],
                          disable_auto_index_creation_and_locking_when_reading_rods=True,
                          max_alternate_alleles=None):

        self.safe_mkdir(splited_dir)

        regions_list,\
            scaffold_to_region_correspondence_dict = self.prepare_region_list_by_length(max_length=max_total_scaffold_length_per_chunk,
                                                                                        max_seq_number=max_scaffold_number_per_chunk,
                                                                                        length_dict=length_dict,
                                                                                        reference=None if length_dict is not None else reference,
                                                                                        parsing_mode=parsing_mode,
                                                                                        output_dir="%s/regions/" % splited_dir,
                                                                                        split_scaffolds=False) if region_list is None else region_list

        options = self.parse_options_for_parallel_run(reference, gvcf_list, extension_list=extension_list,
                                                      max_alternate_alleles=max_alternate_alleles,
                                                      disable_auto_index_creation_and_locking_when_reading_rods=disable_auto_index_creation_and_locking_when_reading_rods)

        output_index = 1
        options_list = []

        region_vcf_list = []

        for regions in regions_list:
            region_options = " -o %s/%s_%i.vcf" % (splited_dir, splited_prefix, output_index)
            region_vcf_list.append("%s/%s_%i.vcf" % (splited_dir, splited_prefix, output_index))
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
                                              output_vcf,
                                              close_fd_after=False,
                                              extension_list=[".vcf", ])

        ValidateVariants.jar_path = self.jar_path
        ValidateVariants.jar = self.jar

        ValidateVariants.index_vcf(reference, output_vcf)
