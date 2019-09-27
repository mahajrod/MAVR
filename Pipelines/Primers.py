#!/usr/bin/env python
from copy import deepcopy
from collections import OrderedDict
from RouToolPa.Tools.Primers import Primer3
from RouToolPa.Tools.Kmers import Glistmaker
from RouToolPa.Tools.RepeatMasking import TRF
from RouToolPa.Parsers.Primer3 import CollectionPrimer3
from RouToolPa.Routines import AnnotationsRoutines
from RouToolPa.Parsers.Sequence import CollectionSequence
from RouToolPa.Parsers.GFF import CollectionGFF
from Pipelines.Abstract import Pipeline


class PrimerPipeline(Pipeline):

    def __init__(self, primer3_dir="", primer3_bin="primer3_core",
                 glistmaker_dir="", glistmaker_bin="glistmaker", primer3_thermo_config_dir=None):
        Pipeline.__init__(self)
        self.primer3_dir = primer3_dir
        self.primer3_bin = primer3_bin
        self.glistmaker_dir = glistmaker_dir
        self.glistmaker_bin = glistmaker_bin
        self.primer3_thermo_config_dir = primer3_thermo_config_dir

    def prepare_primer3_files(*args, **kwargs):
        # Pipeline-specific function
        pass

    def predict_primers(self, coordinates_file, fasta_with_flanks, output_prefix,
                        directory_with_kmer_counts, kmer_file_prefix, pcr_product_size_range=None,
                        optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                        softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                        optimal_melting_temperature=None, min_melting_temperature=None,
                        max_melting_temperature=None, black_list_of_seqs_fasta=None,
                        thermodynamic_parameters_dir=None, format_output=None,
                        coordinates_format="bed",
                        relative_core_seq_coords_relative_entry="core_seq_coords_relative"
                        ):

        primer3_config_file = "%s.primer3.config" % output_prefix
        primer3_input_file = "%s.primer3.input" % output_prefix
        primer3_output_file = "%s.primer3.out" % output_prefix
        primer3_error_file = "%s.primer3.error" % output_prefix

        self.prepare_primer3_files(coordinates_file, fasta_with_flanks, primer3_config_file, primer3_input_file,
                                   self.check_dir_path(directory_with_kmer_counts), kmer_file_prefix,
                                   pcr_product_size_range=pcr_product_size_range,
                                   optimal_primer_len=optimal_primer_len,
                                   min_primer_len=min_primer_len, max_primer_len=max_primer_len,
                                   max_ns_accepted=max_ns_accepted,
                                   softmasked_input=softmasked_input,
                                   optimal_GC=optimal_GC, min_GC=min_GC, max_GC=max_GC,
                                   optimal_melting_temperature=optimal_melting_temperature,
                                   min_melting_temperature=min_melting_temperature,
                                   max_melting_temperature=max_melting_temperature,
                                   black_list_of_seqs_fasta=black_list_of_seqs_fasta,
                                   thermodynamic_parameters_dir=thermodynamic_parameters_dir if thermodynamic_parameters_dir else self.primer3_thermo_config_dir,
                                   coordinates_format=coordinates_format)
        Primer3.path = self.primer3_dir
        Primer3.predict_primers(primer3_input_file,
                                output_file=primer3_output_file,
                                settings_file=primer3_config_file,
                                error_file=primer3_error_file,
                                format_output=format_output, strict_tags=True
                                )


class STRPrimerPipeline(PrimerPipeline):

    def __init__(self, primer3_dir="", trf_dir="", primer3_bin="primer3_core", trf_bin="trf",
                 glistmaker_dir="", glistmaker_bin="glistmaker", primer3_thermo_config_dir=None):

        PrimerPipeline.__init__(self, primer3_dir=primer3_dir,
                                primer3_bin=primer3_bin,
                                glistmaker_dir=glistmaker_dir,
                                glistmaker_bin=glistmaker_bin,
                                primer3_thermo_config_dir=primer3_thermo_config_dir)
        self.trf_dir = trf_dir
        self.trf_bin = trf_bin

    def prepare_primer3_files(self, trf_flank_gff, fasta_with_flanks, primer3_config_file, primer3_input_file,
                              directory_with_kmer_counts, kmer_file_prefix, pcr_product_size_range=None,
                              optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                              softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                              optimal_melting_temperature=None, min_melting_temperature=None,
                              max_melting_temperature=None, black_list_of_seqs_fasta=None,
                              thermodynamic_parameters_dir=None, relative_core_seq_coords_relative_entry="core_seq_coords_relative"
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
                             directory_with_kmer_counts=directory_with_kmer_counts, kmer_file_prefix=kmer_file_prefix,
                             thermodynamic_parameters_dir=thermodynamic_parameters_dir)

        sequence_dict = self.parse_seq_file(fasta_with_flanks, "parse", format="fasta")

        number_of_sequences = len(sequence_dict)

        with open(primer3_input_file, "w") as primer3_in_fd:
            with open(trf_flank_gff, "r") as trf_gff_fd:
                for line in trf_gff_fd:
                    if line[0] == "#":
                        continue
                    description_dict = self.get_description_dict_from_gff_string(line)

                    repeat_id = description_dict["ID"]
                    # convert_coordinates to 0-based
                    coordinates = map(lambda s: int(s) - 1, description_dict[relative_core_seq_coords_relative_entry].split(","))
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
                                                                        target_list=((coordinates[0], repeat_length),),
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
                        thermodynamic_parameters_dir=None, format_output=None,
                        coordinates_format="gff",
                        relative_core_seq_coords_relative_entry="core_seq_coords_relative"):

        primer3_config_file = "%s.primer3.config" % output_prefix
        primer3_input_file = "%s.primer3.input" % output_prefix
        primer3_output_file = "%s.primer3.out" % output_prefix
        primer3_error_file = "%s.primer3.error" % output_prefix

        self.prepare_primer3_files(trf_flank_gff, fasta_with_flanks, primer3_config_file, primer3_input_file,
                                   self.check_dir_path(directory_with_kmer_counts), kmer_file_prefix,
                                   pcr_product_size_range=pcr_product_size_range,
                                   optimal_primer_len=optimal_primer_len,
                                   min_primer_len=min_primer_len, max_primer_len=max_primer_len,
                                   max_ns_accepted=max_ns_accepted,
                                   softmasked_input=softmasked_input,
                                   optimal_GC=optimal_GC, min_GC=min_GC, max_GC=max_GC,
                                   optimal_melting_temperature=optimal_melting_temperature,
                                   min_melting_temperature=min_melting_temperature,
                                   max_melting_temperature=max_melting_temperature,
                                   black_list_of_seqs_fasta=black_list_of_seqs_fasta,
                                   thermodynamic_parameters_dir=thermodynamic_parameters_dir if thermodynamic_parameters_dir else self.primer3_thermo_config_dir,
                                   relative_core_seq_coords_relative_entry=relative_core_seq_coords_relative_entry
                                   )
        Primer3.path = self.primer3_dir
        Primer3.predict_primers(primer3_input_file,
                                output_file=primer3_output_file,
                                settings_file=primer3_config_file,
                                error_file=primer3_error_file,
                                format_output=format_output, strict_tags=True
                                )

    def primer_prediction_pipeline(self, genome_fasta, output_prefix, trf_gff=None, min_str_period=3, max_str_period=5,
                                   min_copy_number=20, max_copy_number=None, pattern=None, min_perfect_copy_number=20,
                                   require_tandem_perfect_copies=True, left_flank_len=200, right_flank_len=200,
                                   core_seq_coords_entry="core_seq_coords", id_description_entry="ID",
                                   kmer_dir=None, kmer_file_prefix=None, count_kmers=False,
                                   min_percentage_of_matches=None, max_percentage_of_indels=None,
                                   optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                                   softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                                   optimal_melting_temperature=None, min_melting_temperature=None,
                                   max_melting_temperature=None, black_list_of_seqs_fasta=None,
                                   trf_matching_weight=2, trf_mismatching_penalty=7,
                                   trf_indel_penalty=7, trf_matching_probability=80, trf_indel_probability=10,
                                   trf_min_score=50, trf_max_period_size=500, threads=None, min_gap_len=5):

        TRF.path = self.trf_dir
        TRF.threads = threads if threads else self.threads

        Primer3.path = self.primer3_dir
        Primer3.threads = threads if threads else self.threads

        Glistmaker.path = self.glistmaker_dir
        Glistmaker.threads = threads if threads else self.threads

        trf_output_gff = "%s.with_rep_seqs.gff" % output_prefix if trf_gff is None else trf_gff

        filtered_suffix = ""
        filtered_suffix += ".min_period_%i" % min_str_period if min_str_period else ""
        filtered_suffix += ".max_period_%i" % max_str_period if max_str_period else ""
        filtered_suffix += ".min_copy_%i" % min_copy_number if min_copy_number else ""
        filtered_suffix += ".max_copy_%i" % max_copy_number if max_copy_number else ""
        filtered_suffix += ".pattern_%s" % pattern if pattern else ""

        filtered_trf_gff = "%s%s.gff" % (output_prefix, filtered_suffix)
        filtered_out_trf_gff = "%s%s.filtered_out.gff" % (output_prefix, filtered_suffix)

        final_filtered_gff = filtered_trf_gff

        if min_perfect_copy_number:
            filtering_prefix = "%s%s.%s" % (output_prefix, filtered_suffix,
                                                            "min_tandem_perfect_copy_%i" % min_perfect_copy_number if require_tandem_perfect_copies else "min_perfect_copy_%i" % min_perfect_copy_number)

            final_filtered_gff = "%s.gff" % filtering_prefix
            filtered_out_exact_copy_trf_gff = "%s.filtered_out.gff" % filtering_prefix
            #final_filtered_gff = filtered_exact_copy_trf_gff

        final_filtered_len_file = "%s.monomer_len.len" % final_filtered_gff[:-4]

        with_flanks_prefix = "%s.with_flanks" % final_filtered_gff[:-4]
        with_flanks_gff = "%s.gff" % with_flanks_prefix
        with_flanks_fasta = "%s.fasta" % with_flanks_prefix

        primer3_output_prefix = "%s.primer3" % with_flanks_prefix

        if trf_gff is None:
            print("Annotating repeats...")
            trf_report = TRF.parallel_search_tandem_repeat(genome_fasta, output_prefix,
                                                           matching_weight=trf_matching_weight,
                                                           mismatching_penalty=trf_mismatching_penalty,
                                                           indel_penalty=trf_indel_penalty,
                                                           match_probability=trf_matching_probability,
                                                           indel_probability=trf_indel_probability,
                                                           min_alignment_score=trf_min_score,
                                                           max_period=trf_max_period_size,
                                                           report_flanking_sequences=False,
                                                           max_len_per_file=1000000,
                                                           store_intermediate_files=False)
        print("Filtering repeats...")
        TRF.filter_trf_gff(trf_output_gff, filtered_trf_gff, filtered_out_trf_gff, min_period=min_str_period,
                           max_period=max_str_period, min_copy_number=min_copy_number, max_copy_number=max_copy_number,
                           pattern=pattern, min_percentage_of_matches=min_percentage_of_matches,
                           max_percentage_of_indels=max_percentage_of_indels, min_entropy=None, max_entropy=None)

        id_based_location_dict = AnnotationsRoutines.get_id_based_dict_from_gff(trf_output_gff,
                                                                                id_entry=id_description_entry) if trf_gff else trf_report.get_id_based_dict()

        if min_perfect_copy_number:
            TRF.filter_trf_gff_by_exact_copy_number(filtered_trf_gff, final_filtered_gff,
                                                    filtered_out_exact_copy_trf_gff, min_perfect_copy_number,
                                                    perfect_tandem=require_tandem_perfect_copies)

        #print final_filtered_gff
        #print final_filtered_len_file
        TRF.get_monomer_len_file_from_trf_gff(final_filtered_gff, final_filtered_len_file)

        monomer_length_id_file_prefix = "%s.monomer_len" % final_filtered_gff[:-4]
        monomer_length_id_dict = self.split_ids_from_len_file_by_len(final_filtered_len_file,
                                                                     monomer_length_id_file_prefix,
                                                                     len_column=1, id_column=0)

        AnnotationsRoutines.add_flanks_to_gff_record(final_filtered_gff, with_flanks_prefix,
                                                     left_flank_len, right_flank_len, genome_fasta,
                                                     coords_description_entry=core_seq_coords_entry,
                                                     id_description_entry=id_description_entry)

        AnnotationsRoutines.extract_sequences_by_gff(genome_fasta,
                                                     with_flanks_gff,
                                                     with_flanks_fasta,
                                                     type_list="repeat",
                                                     parsing_mode="parse",
                                                     format="fasta")

        if count_kmers:
            print("Counting kmers...")
            if (not kmer_file_prefix) or (not kmer_dir):
                raise ValueError("No kmer file prefix of kmer directory was set")
            glistmaker_prefix = "%s/%s" % (kmer_dir, kmer_file_prefix)
            self.safe_mkdir(kmer_dir)
            Glistmaker.generate_kmer_lists_for_primer3(genome_fasta, glistmaker_prefix, threads=None,
                                                       max_tmp_table_number=None, max_tmp_table_size=None)
        print("Generating primers...")
        for human_readable_output in False, True:
            output_file_prefix = "%s.human_readable" % with_flanks_prefix if human_readable_output else with_flanks_prefix
            self.predict_primers(with_flanks_gff, with_flanks_fasta, output_file_prefix,
                                 kmer_dir, kmer_file_prefix, pcr_product_size_range=None,
                                 optimal_primer_len=optimal_primer_len,
                                 min_primer_len=min_primer_len, max_primer_len=max_primer_len,
                                 max_ns_accepted=max_ns_accepted,
                                 softmasked_input=softmasked_input,
                                 optimal_GC=optimal_GC, min_GC=min_GC, max_GC=max_GC,
                                 optimal_melting_temperature=optimal_melting_temperature,
                                 min_melting_temperature=min_melting_temperature,
                                 max_melting_temperature=max_melting_temperature,
                                 black_list_of_seqs_fasta=black_list_of_seqs_fasta,
                                 thermodynamic_parameters_dir=self.primer3_thermo_config_dir,
                                 format_output=human_readable_output,
                                 relative_core_seq_coords_relative_entry="%s_relative" % core_seq_coords_entry)

        primer3_output_file = "%s.out" % primer3_output_prefix

        filtered_results_file = "%s.filtered.res" % primer3_output_prefix
        filtered_results_table_form_file = "%s.filtered.table_form.res" % primer3_output_prefix
        filtered_results_table_form_with_aln_file = "%s.filtered.table_form_with_aln.res" % primer3_output_prefix
        filtered_out_results_file = "%s.filtered_out.res" % primer3_output_prefix

        primer3_results = CollectionPrimer3(primer3_file=primer3_output_file, from_file=True, id_based_location_dict=id_based_location_dict)

        primer3_results.remove_primers_with_gaps_in_pcr_product(min_gap_len)
        primer3_filtered_results, primer3_filtered_out_results = primer3_results.filter_out_records_without_primers()

        primer3_filtered_results.write(filtered_results_file)
        primer3_filtered_results.write_table_form(filtered_results_table_form_file)
        primer3_filtered_results.write_table_form_with_alignments(filtered_results_table_form_with_aln_file)
        primer3_filtered_out_results.write(filtered_out_results_file)

        filtered_results_file_splited_by_len_prefix = "%s.filtered.monomer_len" % primer3_output_prefix

        stat_fd = open("%s.stats" % output_prefix, "w")

        sorted_monomer_length_list = map(str, sorted(map(int, monomer_length_id_dict.keys())))

        for monomer_length in sorted_monomer_length_list:
            primer3_monomer_len_results = primer3_filtered_results.extract_records_by_ids(monomer_length_id_dict[monomer_length])
            primer3_monomer_len_results.write("%s.%s.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))
            primer3_monomer_len_results.write_table_form("%s.%s.table_form.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))
            primer3_monomer_len_results.write_table_form_with_alignments("%s.%s.table_form_with_aln.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))

            primer3_monomer_len_results.write_table_form2("%s.%s.table_form2.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))
            primer3_monomer_len_results.write_table_form2_short("%s.%s.table_form2_short.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))

            stat_string = "STR monomer length %s bp: %i repeats with primers" % (str(monomer_length), len(primer3_monomer_len_results.records))
            print(stat_string)

            stat_fd.write(stat_string + "\n")

        stat_fd.close()


class MitochondrialAmplificationPrimerPipeline(PrimerPipeline):

    def prepare_primer3_files(self, coordinates_file, sequence_file, primer3_config_file, primer3_input_file,
                              directory_with_kmer_counts, kmer_file_prefix, pcr_product_size_range=(4200, 4700),
                              optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=0,
                              softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                              optimal_melting_temperature=69, min_melting_temperature=68,
                              max_melting_temperature=70, black_list_of_seqs_fasta=None,
                              thermodynamic_parameters_dir=None,
                              coordinates_format="bed",
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
                             directory_with_kmer_counts=directory_with_kmer_counts, kmer_file_prefix=kmer_file_prefix,
                             thermodynamic_parameters_dir=thermodynamic_parameters_dir)

        sequences = CollectionSequence(sequence_file, parsing_mode="parse", format="fasta")
        sequences.get_stats_and_features(count_gaps=False, sort=False)
        half_rotated_sequences = deepcopy(sequences)
        rotation_dict = OrderedDict()
        print sequences.seq_lengths
        for record_id in half_rotated_sequences.records:
            rotation_dict[record_id] = int(sequences.seq_lengths["length"].loc[record_id] / 2)
            half_rotated_sequences.records[record_id] = half_rotated_sequences.records[record_id][rotation_dict[record_id]:] + half_rotated_sequences.records[record_id][:rotation_dict[record_id]]

        coordinates = CollectionGFF(coordinates_file, format=coordinates_format)

        record_config = Primer3.generate_config(primer_task=None,
                                                pick_left_primer=True,
                                                pick_right_primer=True,
                                                pick_internal_oligo=False,
                                                optimal_primer_len=None,
                                                min_primer_len=None, max_primer_len=None,
                                                max_ns_accepted=None,
                                                pcr_product_size_range=None,
                                                explain_primers=None, report_all_primers=None,
                                                softmasked_input=False,
                                                optimal_GC=None, min_GC=None, max_GC=None,
                                                optimal_melting_temperature=None,
                                                min_melting_temperature=None,
                                                max_melting_temperature=None,
                                                black_list_of_seqs_fasta=None, mask_sequence=False,
                                                directory_with_kmer_counts=None, kmer_file_prefix=None)

        with open(primer3_input_file, "w") as primer3_in_fd:
            for record_id in sequences.records:
                index = 1
                for line_tuple in coordinates.records.loc[[record_id]].itertuples(index=False):
                    start = getattr(line_tuple, "start")
                    end = getattr(line_tuple, "end")
                    seq = sequences.records[record_id]
                    if start > end:
                        start, end = start - rotation_dict[record_id], end + rotation_dict[record_id]
                        seq = half_rotated_sequences.records[record_id]

                    primer_input_record = Primer3.generate_input_record("%s.primer_%i" % (record_id, index),
                                                                        seq,
                                                                        target_list=((start, end - start),),
                                                                        excluded_region_list=None,
                                                                        included_region_list=None,
                                                                        internal_oligo_excluded_region_list=None,
                                                                        overlap_junction_list=None,
                                                                        forward_primer=None,
                                                                        reverse_primer=None,
                                                                        internal_oligo=None,
                                                                        force_forward_start=None,
                                                                        force_forward_end=None,
                                                                        force_reverse_start=None,
                                                                        force_reverse_end=None,
                                                                        config=record_config)
                    primer3_in_fd.write(primer_input_record)

                    index += 1

    def predict_primers(self, coordinates_gff, fasta_file, output_prefix,
                        directory_with_kmer_counts, kmer_file_prefix, pcr_product_size_range=(4200, 4700),
                        optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                        softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                        optimal_melting_temperature=69, min_melting_temperature=68,
                        max_melting_temperature=70, black_list_of_seqs_fasta=None,
                        thermodynamic_parameters_dir=None, format_output=None,
                        coordinates_format="bed",
                        relative_core_seq_coords_relative_entry="core_seq_coords_relative"):

        primer3_config_file = "%s.primer3.config" % output_prefix
        primer3_input_file = "%s.primer3.input" % output_prefix
        primer3_output_file = "%s.primer3.out" % output_prefix
        primer3_error_file = "%s.primer3.error" % output_prefix

        self.prepare_primer3_files(coordinates_gff, fasta_file, primer3_config_file, primer3_input_file,
                                   self.check_dir_path(directory_with_kmer_counts), kmer_file_prefix,
                                   pcr_product_size_range=pcr_product_size_range,
                                   optimal_primer_len=optimal_primer_len,
                                   min_primer_len=min_primer_len, max_primer_len=max_primer_len,
                                   max_ns_accepted=max_ns_accepted,
                                   softmasked_input=softmasked_input,
                                   optimal_GC=optimal_GC, min_GC=min_GC, max_GC=max_GC,
                                   optimal_melting_temperature=optimal_melting_temperature,
                                   min_melting_temperature=min_melting_temperature,
                                   max_melting_temperature=max_melting_temperature,
                                   black_list_of_seqs_fasta=black_list_of_seqs_fasta,
                                   thermodynamic_parameters_dir=thermodynamic_parameters_dir if thermodynamic_parameters_dir else self.primer3_thermo_config_dir,
                                   )
        Primer3.path = self.primer3_dir
        Primer3.predict_primers(primer3_input_file,
                                output_file=primer3_output_file,
                                settings_file=primer3_config_file,
                                error_file=primer3_error_file,
                                format_output=format_output, strict_tags=True
                                )

    def primer_prediction_pipeline(self, coordinates_gff, mt_fasta, genome_fasta, output_prefix,
                                   kmer_dir=None, kmer_file_prefix=None, pcr_product_size_range=(4200, 4700),
                                   count_kmers=False,
                                   optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                                   softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                                   optimal_melting_temperature=69, min_melting_temperature=68,
                                   max_melting_temperature=70, black_list_of_seqs_fasta=None,
                                   threads=None,):

        Primer3.path = self.primer3_dir
        Primer3.threads = threads if threads else self.threads

        Glistmaker.path = self.glistmaker_dir
        Glistmaker.threads = threads if threads else self.threads
        sequences = CollectionSequence(mt_fasta, parsing_mode="parse", format="fasta")
        if count_kmers:
            print("Counting kmers...")
            if (not kmer_file_prefix) or (not kmer_dir):
                raise ValueError("No kmer file prefix of kmer directory was set")
            glistmaker_prefix = "%s/%s" % (kmer_dir, kmer_file_prefix)
            self.safe_mkdir(kmer_dir)
            Glistmaker.generate_kmer_lists_for_primer3(genome_fasta, glistmaker_prefix, threads=None,
                                                       max_tmp_table_number=None, max_tmp_table_size=None)
        print("Generating primers...")
        for human_readable_output in False, True:
            output_file_prefix = "%s.human_readable" % output_prefix if human_readable_output else output_prefix
            self.predict_primers(coordinates_gff, mt_fasta, output_file_prefix,
                                 kmer_dir, kmer_file_prefix, pcr_product_size_range=pcr_product_size_range,
                                 optimal_primer_len=optimal_primer_len,
                                 min_primer_len=min_primer_len, max_primer_len=max_primer_len,
                                 max_ns_accepted=max_ns_accepted,
                                 softmasked_input=softmasked_input,
                                 optimal_GC=optimal_GC, min_GC=min_GC, max_GC=max_GC,
                                 optimal_melting_temperature=optimal_melting_temperature,
                                 min_melting_temperature=min_melting_temperature,
                                 max_melting_temperature=max_melting_temperature,
                                 black_list_of_seqs_fasta=black_list_of_seqs_fasta,
                                 thermodynamic_parameters_dir=self.primer3_thermo_config_dir,
                                 format_output=human_readable_output)
        primer3_output_file = "%s.primer3.out" % output_prefix
        primer3_results = CollectionPrimer3(primer3_file=primer3_output_file, from_file=True,
                                            id_based_location_dict=None)
        primer3_results.correct_coordinates(sequences.records)
        primer3_results.write_table_form("%s.table_form.res" % output_prefix)
        primer3_results.write_table_form2("%s.table_form2.res" % output_prefix)
        primer3_results.write_table_form2_short("%s.table_form2_short.res" % output_prefix)

        """
        primer3_output_file = "%s.out" % primer3_output_prefix

        filtered_results_file = "%s.filtered.res" % primer3_output_prefix
        filtered_results_table_form_file = "%s.filtered.table_form.res" % primer3_output_prefix
        filtered_results_table_form_with_aln_file = "%s.filtered.table_form_with_aln.res" % primer3_output_prefix
        filtered_out_results_file = "%s.filtered_out.res" % primer3_output_prefix

        primer3_results = CollectionPrimer3(primer3_file=primer3_output_file, from_file=True, id_based_location_dict=id_based_location_dict)

        primer3_results.remove_primers_with_gaps_in_pcr_product(min_gap_len)
        primer3_filtered_results, primer3_filtered_out_results = primer3_results.filter_out_records_without_primers()

        primer3_filtered_results.write(filtered_results_file)
        primer3_filtered_results.write_table_form(filtered_results_table_form_file)
        primer3_filtered_results.write_table_form_with_alignments(filtered_results_table_form_with_aln_file)
        primer3_filtered_out_results.write(filtered_out_results_file)

        filtered_results_file_splited_by_len_prefix = "%s.filtered.monomer_len" % primer3_output_prefix

        stat_fd = open("%s.stats" % output_prefix, "w")

        sorted_monomer_length_list = map(str, sorted(map(int, monomer_length_id_dict.keys())))

        for monomer_length in sorted_monomer_length_list:
            primer3_monomer_len_results = primer3_filtered_results.extract_records_by_ids(monomer_length_id_dict[monomer_length])
            primer3_monomer_len_results.write("%s.%s.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))
            primer3_monomer_len_results.write_table_form("%s.%s.table_form.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))
            primer3_monomer_len_results.write_table_form_with_alignments("%s.%s.table_form_with_aln.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))

            primer3_monomer_len_results.write_table_form2("%s.%s.table_form2.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))
            primer3_monomer_len_results.write_table_form2_short("%s.%s.table_form2_short.res" % (filtered_results_file_splited_by_len_prefix, monomer_length))

            stat_string = "STR monomer length %s bp: %i repeats with primers" % (str(monomer_length), len(primer3_monomer_len_results.records))
            print(stat_string)

            stat_fd.write(stat_string + "\n")

        stat_fd.close()
        """
"""
~/Soft/MAVR/scripts/repeatmasking/tandem_repeat_masking.py -i ../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta -o assembly.hybrid.all -t 30 -p ~/Soft/TRF/trf

~/Soft/MAVR/scripts/repeatmasking/filter_trf_gff.py -i assembly.hybrid.all.trf.with_rep_seqs.gff -o assembly.hybrid.all.trf.with_rep_seqs.monomer_4.copy_no_less_20.gff -x assembly.hybrid.all.trf.with_rep_seqs.monomer_4.copy_no_less_20.filtered_out.gff -n 4 -m 4 -b 20

~/Soft/MAVR/scripts/repeatmasking/filter_trf_gff_by_exact_copy_number.py -i assembly.hybrid.all.trf.with_rep_seqs.monomer_4.copy_no_less_20.gff -o assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_tandem_copy_no_less_20.gff -x assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_tandem_copy_no_less_20.filtered_out.gff -b 20 -p

~/Soft/MAVR/scripts/repeatmasking/filter_trf_gff_by_exact_copy_number.py -i assembly.hybrid.all.trf.with_rep_seqs.monomer_4.copy_no_less_20.gff -o assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_copy_no_less_20.gff -x assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_copy_no_less_20.filtered_out.gff -b 20

~/Soft/MAVR/scripts/annotation/gff/add_flanks_to_gff_record.py  -i assembly.hybrid.all.trf.with_rep_seqs.monomer_4.exact_tandem_copy_no_less_20.gff -o assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks -f ../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta

~/Soft/MAVR/scripts/sequence/extract_sequences_by_gff.py -i ../../../../assemblies/bionano/assemblies/hybrid_assembly/assembly.hybrid.all.fasta -g assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.gff -o assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.fasta -t repeat

~/Soft/MAVR/scripts/primers/predict_str_primers.py -f assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.gff -s assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.fasta -o assembly.hybrid.all.trf -k ../../../../assemblies/bionano/assemblies/hybrid_assembly/kmers/ -r mustela_nigripes_hybrid -p ~/Soft/primer3/src/ -y ~/Soft/primer3/src/primer3_config/ -m

~/Soft/MAVR/scripts/primers/predict_str_primers.py -f assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.gff -s assembly.hybrid.all.trf.exact_tandem_copy_no_less_20.with_flanks.fasta -o assembly.hybrid.all.trf -k ../../../../assemblies/bionano/assemblies/hybrid_assembly/kmers/ -r mustela_nigripes_hybrid -p ~/Soft/primer3/src/ -y ~/Soft/primer3/src/primer3_config/
"""


