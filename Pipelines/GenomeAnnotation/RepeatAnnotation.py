#!/usr/bin/env python
import os
from RouToolPa.Tools.BLAST import Windowmasker
from RouToolPa.Tools.RepeatMasking import TRF, RepeatMasker
from Pipelines.Filtering import FilteringPipeline


class RepeatAnnotation(FilteringPipeline):

    def __init__(self):
        FilteringPipeline.__init__(self)

    def prepare_repeat_directories(self, output_directory, repeatmasker=True, trf=True, windowmasker=True):

        repeatmasker_dir = "%s/repeatmasker/" % output_directory if repeatmasker else None
        windowmasker_dir = "%s/windowmasker/" % output_directory if windowmasker else None
        trf_dir = "%s/trf/" % output_directory if trf else None

        for directory in (output_directory, repeatmasker_dir, windowmasker_dir, trf_dir):
            if directory is not None:
                self.safe_mkdir(directory)

        return output_directory, repeatmasker_dir, windowmasker_dir, trf_dir

    def annotate_repeats(self, input_fasta, output_directory, output_prefix,
                         repeatmasker=True, trf=True, windowmasker=True, threads=1,
                         trf_matching_weight=2, trf_mismatching_penalty=7, trf_indel_penalty=7,
                         trf_matching_probability=80, trf_indel_probability=10, trf_min_score=50,
                         trf_max_period_size=500, trf_max_seq_len=100000, trf_store_intermediate_files=False,
                         trf_binary_path="",
                         repeatmasker_soft_masking=True, repeatmasker_engine=None, repeatmasker_search_speed=None,
                         repeatmasker_no_low_complexity=None, repeatmasker_only_low_complexity=None,
                         repeatmasker_no_interspersed=None, repeatmasker_only_interspersed=None,
                         repeatmasker_no_rna=None, repeatmasker_only_alu=None, repeatmasker_custom_library=None,
                         repeatmasker_species=None, repeatmasker_html_output=False, repeatmasker_ace_output=False,
                         repeatmasker_gff_output=False):

        if "/" in output_prefix:
            raise ValueError("ERROR!!! Presence of '/' in output prefix. "
                             "Output prefix should be only local prefix, without directories.")

        current_dir = os.getcwd()
        masking_dir, repeatmasker_dir, windowmasker_dir, trf_dir = \
            self.prepare_repeat_directories(output_directory, repeatmasker=repeatmasker,
                                            trf=trf, windowmasker=windowmasker)

        trf_prefix = "%s/%s.trf" % (trf_dir, output_prefix)
        windowmasker_prefix = "%s/%s" % (windowmasker_dir, output_prefix)

        Windowmasker.masking(input_fasta, windowmasker_prefix, input_format="fasta", counts_format="obinary",
                             masking_format="interval", source="windowmasker", feature_type="repeat")

        TRF.threads = threads
        RepeatMasker.threads = threads
        if trf_binary_path:
            trf_path_list = self.split_filename(trf_binary_path)

            TRF.path = trf_path_list[0]
            TRF.cmd = trf_path_list[1] + (trf_path_list[2] if trf_path_list[2] else "")

        TRF.parallel_search_tandem_repeat(input_fasta, trf_prefix, matching_weight=trf_matching_weight,
                                          mismatching_penalty=trf_mismatching_penalty,
                                          indel_penalty=trf_indel_penalty,
                                          match_probability=trf_matching_probability,
                                          indel_probability=trf_indel_probability, min_alignment_score=trf_min_score,
                                          max_period=trf_max_period_size,
                                          report_flanking_sequences=False,
                                          max_len_per_file=trf_max_seq_len,
                                          store_intermediate_files=trf_store_intermediate_files)

        repeatmasker_prefix = "%s/%s%s" % (repeatmasker_dir, self.split_filename(input_fasta)[1],
                                           self.split_filename(input_fasta)[2])
        repeatmasker_out_file = "%s.out" % repeatmasker_prefix # self.split_filename(input_fasta)[1] + self.split_filename(input_fasta)[2])

        RepeatMasker.mask(input_fasta, output_dir=repeatmasker_dir, soft_masking=repeatmasker_soft_masking,
                          engine=repeatmasker_engine,
                          search_speed=repeatmasker_search_speed,
                          no_low_complexity=repeatmasker_no_low_complexity,
                          only_low_complexity=repeatmasker_only_low_complexity,
                          no_interspersed=repeatmasker_no_interspersed,
                          only_interspersed=repeatmasker_only_interspersed,
                          no_rna=repeatmasker_no_rna,
                          only_alu=repeatmasker_only_alu,
                          custom_library=repeatmasker_custom_library,
                          species=repeatmasker_species,
                          html_output=repeatmasker_html_output,
                          ace_output=repeatmasker_ace_output,
                          gff_output=repeatmasker_gff_output)

        repeatmasker_converted_prefix = "%s/%s.repeatmasker" % (repeatmasker_dir, output_prefix)

        repeatmasker_converted_gff = "%s.gff" % repeatmasker_converted_prefix
        repeatmasker_repeat_classes_file = "%s.repeat_classes" % repeatmasker_converted_prefix
        repeatmasker_converted_repeat_families_file = "%s.repeat_families" % repeatmasker_converted_prefix

        RepeatMasker.convert_rm_out_to_gff(repeatmasker_out_file,
                                           repeatmasker_converted_gff,
                                           repeatmasker_repeat_classes_file,
                                           repeatmasker_converted_repeat_families_file)

        merged_output = "%s/%s.repeatmasker.trf.windowmasker.gff" % (output_directory, output_prefix)
        merge_cmd = "sort -k1,1 -k4,4n -k5,5 %s.gff %s.gff %s.windowmasker.gff > %s" % (repeatmasker_converted_prefix,
                                                                           trf_prefix,
                                                                           windowmasker_prefix,
                                                                           merged_output)

        os.system(merge_cmd)
