#!/usr/bin/env python

import os
import shutil

from collections import OrderedDict

from Tools.Kmers import Glistmaker
from Tools.Primers import Primer3
from Tools.RepeatMasking import TRF

from Pipelines.Abstract import Pipeline


class STRPrimerPipeline(Pipeline):

    def __init__(self):
        Pipeline.__init__(self)

    def prepare_primer3_files(self, trf_flank_gff, fasta_with_flanks, primer3_config_file, primer3_input_file,
                              directory_with_kmer_counts, kmer_file_prefix, pcr_product_size_range=None,
                              optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                              softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                              optimal_melting_temperature=None, min_melting_temperature=None,
                              max_melting_temperature=None, black_list_of_seqs_fasta=None,
                              ):

        Primer3.write_config(primer3_config_file, primer_task="generic", pick_left_primer=True, pick_right_primer=True,
                             pick_internal_oligo=False, optimal_primer_len=optimal_primer_len,
                             min_primer_len=min_primer_len, max_primer_len=max_primer_len,
                             max_ns_accepted=max_ns_accepted, pcr_product_size_range=pcr_product_size_range,
                             explain_primers=True, report_all_primers=True,
                             softmasked_input=softmasked_input, optimal_GC=optimal_GC, min_GC=min_GC, max_GC=max_GC,
                             optimal_melting_temperature=optimal_melting_temperature,
                             min_melting_temperature=min_melting_temperature,
                             max_melting_temperature=max_melting_temperature,
                             black_list_of_seqs_fasta=black_list_of_seqs_fasta, mask_sequence=True,
                             directory_with_kmer_counts=directory_with_kmer_counts, kmer_file_prefix=kmer_file_prefix)

        sequence_dict = self.parse_seq_file(fasta_with_flanks, "parse", format="fasta")

        number_of_sequences = len(sequence_dict)

        with open(primer3_input_file, "w") as primer3_in_fd:
            with open(trf_flank_gff, "r") as trf_gff_fd:
                for line in trf_gff_fd:
                    if line[0] == "#":
                        continue
                    description_dict = self.get_description_dict_from_gff_string(line)

                    repeat_id = description_dict["ID"]
                    coordinates = map(int, description_dict["core_seq_coords_relative"].split(","))
                    repeat_length = coordinates[1] - coordinates[0] + 1

                    pcr_product_min_size = int(repeat_length + 20)
                    pcr_product_max_size = pcr_product_min_size * 2

                    record_config = Primer3.generate_config(primer_task=None,
                                                            pick_left_primer=True,
                                                            pick_right_primer=True,
                                                            pick_internal_oligo=False,
                                                            optimal_primer_len=None,
                                                            min_primer_len=None, max_primer_len=None,
                                                            max_ns_accepted=None,
                                                            pcr_product_size_range=(pcr_product_min_size, pcr_product_max_size),
                                                            explain_primers=None, report_all_primers=None,
                                                            softmasked_input=False,
                                                            optimal_GC=None, min_GC=None, max_GC=None,
                                                            optimal_melting_temperature=None,
                                                            min_melting_temperature=None,
                                                            max_melting_temperature=None,
                                                            black_list_of_seqs_fasta=None, mask_sequence=False,
                                                            directory_with_kmer_counts=None, kmer_file_prefix=None)

                    primer_input_record = Primer3.generate_input_record(repeat_id, sequence_dict[repeat_id].seq,
                                                                        target_list=coordinates,
                                                                        excluded_region_list=None,
                                                                        included_region_list=None,
                                                                        internal_oligo_excluded_region_list=None,
                                                                        overlap_junction_list=None,
                                                                        forward_primer=None,
                                                                        reverse_primer=None,
                                                                        internal_oligo=None,
                                                                        force_forward_start=None, force_forward_end=None,
                                                                        force_reverse_start=None, force_reverse_end=None,
                                                                        config=record_config)
                    primer3_in_fd.write(primer_input_record)

    def predict_primers(self, trf_flank_gff, fasta_with_flanks, output_prefix,
                              directory_with_kmer_counts, kmer_file_prefix, pcr_product_size_range=None,
                              optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                              softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                              optimal_melting_temperature=None, min_melting_temperature=None,
                              max_melting_temperature=None, black_list_of_seqs_fasta=None,
                              ):

        primer3_config_file = "%s.primer3.config" % output_prefix
        primer3_input_file = "%s.primer3.input" % output_prefix
        primer3_output_file = "%s.primer3.out" % output_prefix
        primer3_error_file = "%s.primer3.error" % output_prefix

        self.prepare_primer3_files(trf_flank_gff, fasta_with_flanks, primer3_config_file, primer3_input_file,
                                   directory_with_kmer_counts, kmer_file_prefix,
                                   pcr_product_size_range=pcr_product_size_range,
                                   optimal_primer_len=optimal_primer_len,
                                   min_primer_len=min_primer_len, max_primer_len=max_primer_len,
                                   max_ns_accepted=max_ns_accepted,
                                   softmasked_input=softmasked_input,
                                   optimal_GC=optimal_GC, min_GC=min_GC, max_GC=max_GC,
                                   optimal_melting_temperature=optimal_melting_temperature,
                                   min_melting_temperature=min_melting_temperature,
                                   max_melting_temperature=max_melting_temperature,
                                   black_list_of_seqs_fasta=black_list_of_seqs_fasta
                                   )

        Primer3.predict_primers(primer3_input_file,
                                output_file=primer3_output_file,
                                settings_file=primer3_config_file,
                                error_file=primer3_error_file,
                                format_output=None, strict_tags=True
                                )




