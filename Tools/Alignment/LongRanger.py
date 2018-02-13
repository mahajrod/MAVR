#!/usr/bin/env python
import os
import shutil
from Tools.Abstract import Tool


class LongRanger(Tool):

    def __init__(self, path="", max_threads=4, max_memory=None):
        Tool.__init__(self, "longranger", path=path, max_threads=max_threads, max_memory=max_memory)

    def parse_wgs_options(self, reference, run_id, fastq_dir=None, sample_prefix=None, description=None,
                          library_name=None, lane_list=None, indice_list=None, project_name=None, variant_calling_mode=None,
                          gatk_jar_path=None, use_somatic_sv_caller=None, precalled_vcf=None, sample_sex=None,
                          variant_calling_only=None, threads=None, max_memory=None):

        options = " wgs"
        options += " --id=%s" % run_id
        if fastq_dir:
            options += " --fastqs=%s" % fastq_dir
        elif sample_prefix:
            options += " --sample=%s" % sample_prefix
        else:
            raise ValueError("Neither directory with fastq files nor sample prefix was set")

        options += " --reference=%s" % reference
        options += " --description=%s" % description if description else ""
        options += " --library=%s" % library_name if library_name else ""

        options += " --lanes=%s" % (lane_list if isinstance(lane_list, str) else ",".join(lane_list)) if lane_list else ""
        options += " --indices=%s" % (indice_list if isinstance(indice_list, str) else ",".join(indice_list)) if indice_list else ""
        options += " --project=%s" % project_name if project_name else ""

        if variant_calling_mode == "gatk":
            options += " --vcmode=gatk:%s" % gatk_jar_path
        else:
            options += " --vcmode=%s" % variant_calling_mode

        options += " --somatic" if use_somatic_sv_caller else ""
        options += " --precalled=%s" % precalled_vcf if precalled_vcf else ""
        options += " --sex=%s" if sample_sex else ""

        options += " --vconly" if variant_calling_only else ""
        options += " --localcores=%i" % (threads if threads else self.threads)
        options += " --localmem=%i" % (str(max_memory) if max_memory else str(self.max_memory)) if max_memory or self.max_memory else ""

        return options

    def run_wgs_analysis(self, reference, run_id, fastq_dir=None, sample_prefix=None, description=None,
                         library_name=None, lane_list=None, indice_list=None, project_name=None, variant_calling_mode=None,
                         gatk_jar_path=None, use_somatic_sv_caller=None, precalled_vcf=None, sample_sex=None,
                         variant_calling_only=None, threads=None, max_memory=None):

        options = self.parse_wgs_options(reference, run_id, fastq_dir=fastq_dir, sample_prefix=sample_prefix,
                                         description=description, library_name=library_name, lane_list=lane_list,
                                         indice_list=indice_list, project_name=project_name,
                                         variant_calling_mode=variant_calling_mode,
                                         gatk_jar_path=gatk_jar_path, use_somatic_sv_caller=use_somatic_sv_caller,
                                         precalled_vcf=precalled_vcf, sample_sex=sample_sex,
                                         variant_calling_only=variant_calling_only, threads=threads,
                                         max_memory=max_memory)

        self.execute(options=options)

    def run_wgs_analysis_from_fastq(self, reference, run_id, fastq_dir, description=None, library_name=None,
                                    lane_list=None, indice_list=None, project_name=None, variant_calling_mode=None,
                                    gatk_jar_path=None, use_somatic_sv_caller=None, precalled_vcf=None, sample_sex=None,
                                    variant_calling_only=None, threads=None, max_memory=None):

        options = self.parse_wgs_options(reference, run_id, fastq_dir=fastq_dir, sample_prefix=None,
                                         description=description, library_name=library_name, lane_list=lane_list,
                                         indice_list=indice_list, project_name=project_name,
                                         variant_calling_mode=variant_calling_mode,
                                         gatk_jar_path=gatk_jar_path, use_somatic_sv_caller=use_somatic_sv_caller,
                                         precalled_vcf=precalled_vcf, sample_sex=sample_sex,
                                         variant_calling_only=variant_calling_only, threads=threads,
                                         max_memory=max_memory)

        self.execute(options=options)
