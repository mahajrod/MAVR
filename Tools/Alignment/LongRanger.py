#!/usr/bin/env python
import os
import shutil
from collections import OrderedDict
from Tools.Abstract import Tool
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

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
        options += " --localmem=%s" % (str(max_memory) if max_memory else str(self.max_memory)) if max_memory or self.max_memory else ""

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

    def prepare_record_dict_reference(self, record_dict, max_scaffold_length=527000000, max_scaffold_number=500,
                                      polyN_len=500):
        length_dict = self.get_lengths(record_dict)
        id_list = length_dict.keys()
        number_of_scaffolds = len(length_dict)
        polyN_insersion = "N" * polyN_len

        output_dict = OrderedDict()
        if number_of_scaffolds <= max_scaffold_number:
            return record_dict

        total_length_of_short_scaffolds_with_insertion = 0
        for i in range(499, number_of_scaffolds):
            total_length_of_short_scaffolds_with_insertion += polyN_len + len(length_dict[id_list[i]])

        if total_length_of_short_scaffolds_with_insertion < max_scaffold_length:
            merged_record = SeqRecord(id="merged_samll_scaffolds", description="merged_records:%i-%i" % (500, number_of_scaffolds),
                                      seq=polyN_insersion.join(map(lambda record_index: record_dict[id_list[record_index]], range(499, number_of_scaffolds))))
            for i in range(0, 499):
                output_dict[id_list[i]] = record_dict[id_list[i]]
            output_dict[merged_record.id] = merged_record

            return output_dict
        
        raise ValueError("Following code have not been completed yet!!!!!!!!!")
        # TODO: following code have not been completed yet!!!!!!!!!
        number_of_merged_records = 0
        total_number_of_records = number_of_scaffolds
        records_to_merge_list = []

        while total_number_of_records > max_scaffold_number:
            tmp_merge_list = []
            tmp_length = 0
            for i in range(number_of_scaffolds, -1, -1):
                if tmp_length + len(length_dict[id_list[i]]) < max_scaffold_length:
                    pass

    def prepare_reference(self, reference, prepared_reference, max_scaffold_length=527000000, max_scaffold_number=500,
                          polyN_len=500):
        """
        Diploid genome — phasing algorithm currently assumes 2 haplotypes
        500 contigs or fewer — if your assembly has more than 500 contigs, concatenate smaller contigs together with 500 N's separating each original contig, until there are fewer than 500 contigs total
        All contigs must be no more than 2^29-1 bp, or 528Mb, in length; this is a limitation of BAM index file format
        All contigs must have no colons or spaces in their names.
        """

        record_dict = self.parse_seq_file(reference, mode="parse")

        prepared_record_dict = self.prepare_record_dict_reference(record_dict,
                                                                  max_scaffold_length=max_scaffold_length,
                                                                  max_scaffold_number=max_scaffold_number,
                                                                  polyN_len=polyN_len)
        SeqIO.write(list(prepared_record_dict), prepared_reference, format='fasta')


