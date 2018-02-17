#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os


from Tools.Abstract import Tool


class Primer3(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "primer3_core", path=path, max_threads=max_threads)

    @staticmethod
    def parse_options(input_file, output_file=None, settings_file=None, format_output=None,
                      strict_tags=True, error_file=None):

        options = " --format_output" if format_output else ""
        options += " --strict_tags" if strict_tags else ""
        options += " --p3_settings_file=%s" % settings_file if settings_file else ""
        options += " --error=%s" % error_file if error_file else ""
        options += " --output=%s" % output_file if output_file else ''

        options += " %s" % input_file

        return options

    @staticmethod
    def generate_config(primer_task="generic", pick_left_primer=True, pick_right_primer=True, pick_internal_oligo=False,
                        optimal_primer_len=None, min_primer_len=None, max_primer_len=None, max_ns_accepted=None,
                        pcr_product_size_range=None, explain_primers=True, report_all_primers=True,
                        softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                        optimal_melting_temperature=None, min_melting_temperature=None,
                        max_melting_temperature=None, black_list_of_seqs_fasta=None, mask_sequence=False,
                        directory_with_kmer_counts=None, kmer_file_prefix=None, thermodynamic_parameters_dir=None):
        """
        PRIMER_TASK:  generic,
                      check_primers,
                      pick_primer_list,
                      pick_cloning_primers,
                      pick_discriminative_primers

        TODO: A lot of PRIMER options were not implemented!!!
        """
        config = ""
        config += "PRIMER_TASK=%s\n" % primer_task if primer_task else ""
        config += "PRIMER_PICK_LEFT_PRIMER=%i\n" % (1 if pick_left_primer else 0)
        config += "PRIMER_PICK_INTERNAL_OLIGO=%i\n" % (1 if pick_internal_oligo else 0)
        config += "PRIMER_PICK_RIGHT_PRIMER=%i\n" % (1 if pick_right_primer else 0)

        config += "PRIMER_OPT_SIZE=%i\n" % optimal_primer_len if optimal_primer_len else ""
        config += "PRIMER_MIN_SIZE=%i\n" % min_primer_len if min_primer_len else ""
        config += "PRIMER_MAX_SIZE=%i\n" % max_primer_len if max_primer_len else ""

        config += "PRIMER_MAX_NS_ACCEPTED=%i\n" % max_ns_accepted if max_ns_accepted else ""
        config += "PRIMER_PRODUCT_SIZE_RANGE=%s\n" % ("-".join(map(str, pcr_product_size_range))) if pcr_product_size_range else ""

        config += "PRIMER_EXPLAIN_FLAG=1\n" if explain_primers else ""
        config += "P3_FILE_FLAG=1\n" if report_all_primers else ""
        config += "PRIMER_LOWERCASE_MASKING=1\n" if softmasked_input else ""

        config += "PRIMER_OPT_GC_PERCENT=%f\n" % optimal_GC if optimal_GC else ""
        config += "PRIMER_MIN_GC=%f\n" % min_GC if min_GC else ""
        config += "PRIMER_MAX_GC=%f\n" % max_GC if max_GC else ""

        config += "PRIMER_OPT_TM=%f\n" % optimal_melting_temperature if optimal_melting_temperature else ""
        config += "PRIMER_MIN_TM=%f\n" % min_melting_temperature if min_melting_temperature else ""
        config += "PRIMER_MAX_TM=%f\n" % max_melting_temperature if max_melting_temperature else ""

        config += "PRIMER_MISPRIMING_LIBRARY=%s\n" % black_list_of_seqs_fasta if black_list_of_seqs_fasta else ""
        config += "PRIMER_MASK_TEMPLATE=1\n" if mask_sequence else ""

        config += "PRIMER_MASK_KMERLIST_PATH=%s\n" % directory_with_kmer_counts if directory_with_kmer_counts else ""
        config += "PRIMER_MASK_KMERLIST_PREFIX=%s\n" % kmer_file_prefix if kmer_file_prefix else ""
        config += "PRIMER_THERMODYNAMIC_PARAMETERS_PATH=%s\n" % thermodynamic_parameters_dir if thermodynamic_parameters_dir else ""

        return config

    def write_config(self, config_file, primer_task="generic", pick_left_primer=True, pick_right_primer=True,
                     pick_internal_oligo=False, optimal_primer_len=None, min_primer_len=None, max_primer_len=None,
                     max_ns_accepted=None, pcr_product_size_range=None, explain_primers=True, report_all_primers=True,
                     softmasked_input=False, optimal_GC=None, min_GC=None, max_GC=None,
                     optimal_melting_temperature=None, min_melting_temperature=None,
                     max_melting_temperature=None, black_list_of_seqs_fasta=None, mask_sequence=False,
                     directory_with_kmer_counts=None, kmer_file_prefix=None,
                     thermodynamic_parameters_dir=None):

        config = "Primer3 File - http://primer3.sourceforge.net\n"
        config += "P3_FILE_TYPE=settings\n"
        config += "\n"
        config += "P3_FILE_ID=STR primers\n"

        config += self.generate_config(primer_task=primer_task, pick_left_primer=pick_left_primer,
                                       pick_right_primer=pick_right_primer, pick_internal_oligo=pick_internal_oligo,
                                       optimal_primer_len=optimal_primer_len, min_primer_len=min_primer_len,
                                       max_primer_len=max_primer_len, max_ns_accepted=max_ns_accepted,
                                       pcr_product_size_range=pcr_product_size_range, explain_primers=explain_primers,
                                       report_all_primers=report_all_primers, softmasked_input=softmasked_input,
                                       optimal_GC=optimal_GC, min_GC=min_GC, max_GC=max_GC,
                                       optimal_melting_temperature=optimal_melting_temperature,
                                       min_melting_temperature=min_melting_temperature,
                                       max_melting_temperature=max_melting_temperature,
                                       black_list_of_seqs_fasta=black_list_of_seqs_fasta,
                                       mask_sequence=mask_sequence,
                                       directory_with_kmer_counts=directory_with_kmer_counts,
                                       kmer_file_prefix=kmer_file_prefix,
                                       thermodynamic_parameters_dir=thermodynamic_parameters_dir)
        config += "=\n"

        with open(config_file, "w") as config_fd:
            config_fd.write(config)

    @staticmethod
    def generate_input_record(sequence_id, sequence,
                              target_list=None, excluded_region_list=None, included_region_list=None,
                              internal_oligo_excluded_region_list=None, overlap_junction_list=None,
                              forward_primer=None, reverse_primer=None, internal_oligo=None,
                              force_forward_start=None, force_forward_end=None,
                              force_reverse_start=None, force_reverse_end=None,
                              config=""):
        """
        TODO:
            NOT IMPLEMENTED
                    SEQUENCE_PRIMER_PAIR_OK_REGION_LIST (semicolon separated list of integer "quadruples"; default empty)

                        This tag allows detailed specification of possible locations of left and right primers in primer pairs.

                        The associated value must be a semicolon-separated list of

                        <left_start>,<left_length>,<right_start>,<right_length>
                        quadruples. The left primer must be in the region specified by <left_start>,<left_length> and the right primer must be in the region specified by <right_start>,<right_length>. <left_start> and <left_length> specify the location of the left primer in terms of the index of the first base in the region and the length of the region. <right_start> and <right_length> specify the location of the right primer in analogous fashion. As seen in the example below, if no integers are specified for a region then the location of the corresponding primer is not constrained.

                        Example:

                            SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=100,50,300,50 ; 900,60,, ; ,,930,100

                    SEQUENCE_START_CODON_POSITION (int; default -2000000)

                        This parameter should be considered EXPERIMENTAL at this point. Please check the output carefully; some erroneous inputs might cause an error in Primer3.

                        Index of the first base of a start codon. This parameter allows Primer3 to select primer pairs to create in-frame amplicons e.g. to create a template for a fusion protein. Primer3 will attempt to select an in-frame left primer, ideally starting at or to the left of the start codon, or to the right if necessary. Negative values of this parameter are legal if the actual start codon is to the left of available sequence. If this parameter is non-negative Primer3 signals an error if the codon at the position specified by this parameter is not an ATG. A value less than or equal to -10^6 indicates that Primer3 should ignore this parameter.

                        Primer3 selects the position of the right primer by scanning right from the left primer for a stop codon. Ideally the right primer will end at or after the stop codon.

                    SEQUENCE_QUALITY (space separated integers; default empty)

                        A list of space separated integers. There must be exactly one integer for each base in SEQUENCE_TEMPLATE if this argument is non-empty. For example, for the sequence ANNTTCA... SEQUENCE_QUALITY might be 45 10 0 50 30 34 50 67 .... High numbers indicate high confidence in the base called at that position and low numbers indicate low confidence in the base call at that position. This parameter is only relevant if you are using a base calling program that provides quality information (for example phred).
        """

        input_record = "SEQUENCE_ID=%s\n" % sequence_id
        input_record += "SEQUENCE_TEMPLATE=%s\n" % str(sequence)
        input_record += "SEQUENCE_TARGET=%s\n" % (" ".join(map(lambda s: "%s,%s" % (str(s[0]), str(s[1])), target_list))) if target_list else ""
        input_record += "SEQUENCE_INCLUDED_REGION=%s\n" % (" ".join(map(lambda s: "%s,%s" % (str(s[0]), str(s[1])), included_region_list))) if included_region_list else ""
        input_record += "SEQUENCE_EXCLUDED_REGION=%s\n" % (" ".join(map(lambda s: "%s,%s" % (str(s[0]), str(s[1])), excluded_region_list))) if excluded_region_list else ""
        input_record += "SEQUENCE_INTERNAL_EXCLUDED_REGION=%s\n" % (" ".join(map(lambda s: "%s,%s" % (str(s[0]), str(s[1])), internal_oligo_excluded_region_list))) if internal_oligo_excluded_region_list else ""
        input_record += "SEQUENCE_OVERLAP_JUNCTION_LIST=%s\n" % (" ".join(map(str, overlap_junction_list))) if overlap_junction_list else ""

        input_record += "SEQUENCE_PRIMER=%s\n" % forward_primer if forward_primer else ""
        input_record += "SEQUENCE_PRIMER_REVCOMP=%s\n" % reverse_primer if reverse_primer else ""
        input_record += "SEQUENCE_INTERNAL_OLIGO=%s\n" % internal_oligo if internal_oligo else ""

        input_record += "SEQUENCE_FORCE_LEFT_START=%i\n" % force_forward_start if force_forward_start else ""
        input_record += "SEQUENCE_FORCE_LEFT_END=%i\n" % force_forward_end if force_forward_end else ""
        input_record += "SEQUENCE_FORCE_RIGHT_START=%i\n" % force_reverse_start if force_reverse_start else ""
        input_record += "SEQUENCE_FORCE_RIGHT_END=%i\n" % force_reverse_end if force_reverse_end else ""

        input_record += config
        input_record += "=\n"

        return input_record

    def predict_primers(self, input_file, output_file=None, settings_file=None, format_output=None, strict_tags=True,
                        error_file=None):
        options = self.parse_options(input_file,
                                     output_file=output_file,
                                     settings_file=settings_file,
                                     format_output=format_output,
                                     strict_tags=strict_tags,
                                     error_file=error_file)

        self.execute(options=options)



    """
    def search_tandem_repeats(self, query_file, matching_weight=2, mismatching_penalty=7, indel_penalty=7,
                              match_probability=80, indel_probability=10, min_alignment_score=50, max_period=500,
                              report_flanking_sequences=False, make_dat_file=True, disable_html_output=True):

        options = " %s" % query_file
        options += self.parse_common_options(matching_weight=matching_weight, mismatching_penalty=mismatching_penalty,
                                             indel_penalty=indel_penalty, match_probability=match_probability,
                                             indel_probability=indel_probability, min_alignment_score=min_alignment_score,
                                             max_period=max_period, report_flanking_sequences=report_flanking_sequences,
                                             make_dat_file=make_dat_file)

        options += " -h" if disable_html_output else ""

        self.execute(options)



    def parallel_search_tandem_repeat(self, query_file, output_prefix, matching_weight=2, mismatching_penalty=7,
                                      indel_penalty=7,
                                      match_probability=80, indel_probability=10, min_alignment_score=50, max_period=500,
                                      report_flanking_sequences=False, splited_fasta_dir="splited_fasta_dir",
                                      splited_result_dir="splited_output", converted_output_dir="converted_output",
                                      max_len_per_file=100000, store_intermediate_files=False):
        work_dir = os.getcwd()
        splited_filename = FileRoutines.split_filename(query_file)
        self.split_fasta_by_seq_len(query_file, splited_fasta_dir, max_len_per_file=max_len_per_file,
                                    output_prefix=splited_filename[1])

        common_options = self.parse_common_options(matching_weight=matching_weight,
                                                   mismatching_penalty=mismatching_penalty,
                                                   indel_penalty=indel_penalty, match_probability=match_probability,
                                                   indel_probability=indel_probability,
                                                   min_alignment_score=min_alignment_score,
                                                   max_period=max_period,
                                                   report_flanking_sequences=report_flanking_sequences,
                                                   make_dat_file=True)
        common_options += " -h"  # suppress html output
        options_list = []
        splited_files = os.listdir(splited_fasta_dir)

        FileRoutines.safe_mkdir(splited_result_dir)
        FileRoutines.safe_mkdir(converted_output_dir)
        os.chdir(splited_result_dir)

        input_dir = splited_fasta_dir if (splited_fasta_dir[0] == "/") or (splited_fasta_dir[0] == "~") \
                    else "../%s" % splited_fasta_dir

        for filename in splited_files:
            file_options = "%s/%s" % (input_dir, filename)
            file_options += common_options
            options_list.append(file_options)

        self.parallel_execute(options_list)
    """


if __name__ == "__main__":
    pass
