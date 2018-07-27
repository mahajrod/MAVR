#!/usr/bin/env python
import os
import shutil
from collections import OrderedDict
from Tools.Abstract import Tool
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq

from CustomCollections.GeneralCollections import SynDict


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
        options += " --sex=%s" % sample_sex if sample_sex else ""

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

    def prepare_record_dict_reference(self, record_dict, coord_file=None, max_scaffold_length=527000000,
                                      max_scaffold_number=500, polyN_len=500, symbols_to_replace_list=[":"]):
        length_dict = self.get_lengths(record_dict)
        id_list = length_dict.keys()
        number_of_scaffolds = len(length_dict)
        polyN_insersion = "N" * polyN_len
        renamed_scaffolds_syn_dict = SynDict()

        #print id_list
        output_dict = OrderedDict()
        if number_of_scaffolds <= max_scaffold_number:
            return record_dict

        total_length_of_short_scaffolds_with_insertion = 0
        for i in range(499, number_of_scaffolds):
            #print length_dict[id_list[i]]
            total_length_of_short_scaffolds_with_insertion += polyN_len + length_dict[id_list[i]]

        if total_length_of_short_scaffolds_with_insertion < max_scaffold_length:
            merged_record = SeqRecord(id="merged_small_scaffolds", description="merged_records:%i-%i" % (500, number_of_scaffolds),
                                      seq=Seq(polyN_insersion.join(map(lambda record_index: str(record_dict[id_list[record_index]].seq),
                                                                   range(499, number_of_scaffolds)))))
            if coord_file:
                with open(coord_file, "w") as coord_fd:
                    coord_fd.write("#id\tstart\tstop\n")
                    coord_fd.write("%s\t1\t%i\n" % (id_list[499], length_dict[id_list[499]]))
                    prev_end = length_dict[id_list[499]]
                    if number_of_scaffolds > 500:
                        for i in range(500, number_of_scaffolds):
                            coord_fd.write("%s\t%i\t%i\n" % (id_list[i],
                                                              prev_end + 1 + polyN_len,
                                                              prev_end + length_dict[id_list[i]] + polyN_len))
                            prev_end += length_dict[id_list[i]] + polyN_len

            #print merged_record
            for i in range(0, 499):
                if ":" not in id_list[i]:
                    output_dict[id_list[i]] = record_dict[id_list[i]]
                else:
                    new_id = id_list[i]
                    for symbol in symbols_to_replace_list:
                        new_id = new_id.replace(symbol, "_")
                    renamed_scaffolds_syn_dict[id_list[i]] = new_id
                    output_dict[new_id] = SeqRecord(seq=record_dict[id_list[i]].seq, id=new_id, description=record_dict[id_list[i]].description)
            output_dict[merged_record.id] = merged_record
            #print output_dict
            return output_dict, renamed_scaffolds_syn_dict

        # TODO: following code have not been completed yet!!!!!!!!!
        raise ValueError("Following code have not been completed yet!!!!!!!!!")

        number_of_merged_records = 0
        total_number_of_records = number_of_scaffolds
        records_to_merge_list = []

        while total_number_of_records > max_scaffold_number:
            tmp_merge_list = []
            tmp_length = 0
            for i in range(number_of_scaffolds, -1, -1):
                if tmp_length + len(length_dict[id_list[i]]) < max_scaffold_length:
                    pass

    def prepare_reference(self, reference, prepared_reference, renamed_scaffolds_syn_file,
                          max_scaffold_length=527000000, max_scaffold_number=500,
                          polyN_len=500, coord_file=None, symbols_to_replace_list=[":"]):
        """
        Diploid genome - phasing algorithm currently assumes 2 haplotypes
        500 contigs or fewer - if your assembly has more than 500 contigs, concatenate smaller contigs together
        with 500 Ns separating each original contig, until there are fewer than 500 contigs total
        All contigs must be no more than 2^29-1 bp, or 528Mb, in length; this is a limitation of BAM index file format
        All contigs must have no colons or spaces in their names.
        """

        record_dict = self.parse_seq_file(reference, mode="parse")

        prepared_record_dict, renamed_scaffolds_syn_dict = self.prepare_record_dict_reference(record_dict,
                                                                                              max_scaffold_length=max_scaffold_length,
                                                                                              max_scaffold_number=max_scaffold_number,
                                                                                              polyN_len=polyN_len,
                                                                                              coord_file=coord_file,
                                                                                              symbols_to_replace_list=symbols_to_replace_list)
        SeqIO.write(self.record_from_dict_generator(prepared_record_dict), prepared_reference, format='fasta')
        renamed_scaffolds_syn_dict.write(renamed_scaffolds_syn_file)
        self.index_reference(prepared_reference)

    def index_reference(self, reference):

        options = " mkref %s" % reference

        self.execute(options=options)
