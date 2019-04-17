#!/usr/bin/env python

import os
from collections import OrderedDict
from RouToolPa.Tools.Alignment import LongRanger
from Pipelines.Abstract import Pipeline


class TenXAlignmentPipeline(Pipeline):

    def __init__(self, max_threads=1, max_memory=10, longranger_dir=""):
        Pipeline.__init__(self, max_threads=max_threads, max_memory=max_memory)
        self.longranger_dir = longranger_dir if not (longranger_dir is None) else ""

    def prepare_directories(self, output_directory, sample_list):
        print(output_directory)
        directories_dict = {output_directory: OrderedDict()}
        for sample in sample_list:
            directories_dict[output_directory][sample] = OrderedDict()

        self.recursive_mkdir(directories_dict)

    def align_and_call(self, samples_directory, output_directory, reference, samples_to_handle=None,
                       variant_calling_mode=None, gatk_jar_path=None, use_somatic_sv_caller=None,
                       precalled_vcf=None, sample_sex=None,
                       variant_calling_only=None, threads=None, max_memory=None, longranger_dir=None):

        start_dir = os.getcwd()
        out_dir_abs = os.path.abspath(output_directory)
        LongRanger.path = longranger_dir if longranger_dir else self.longranger_dir
        LongRanger.threads = threads if threads else self.threads
        LongRanger.max_memory = max_memory if max_memory else self.max_memory

        sample_list = samples_to_handle if samples_to_handle else self.get_sample_list(samples_directory)

        self.prepare_directories(out_dir_abs, sample_list)

        os.chdir(out_dir_abs)
        for sample in sample_list:
            fastq_dir = "%s/%s/" % (samples_directory, sample)

            LongRanger.run_wgs_analysis(reference, sample, fastq_dir, description=sample, library_name=sample,
                                        lane_list=None, indice_list=None, project_name=None,
                                        variant_calling_mode=variant_calling_mode,
                                        gatk_jar_path=gatk_jar_path, use_somatic_sv_caller=use_somatic_sv_caller,
                                        precalled_vcf=precalled_vcf, sample_sex=sample_sex,
                                        variant_calling_only=variant_calling_only,
                                        max_memory=max_memory if max_memory else self.max_memory)

        os.chdir(start_dir)
