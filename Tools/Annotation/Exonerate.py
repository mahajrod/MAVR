#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os
import shutil

from Tools.Abstract import Tool

from Parsers.TRF import CollectionTRF
from Routines.File import split_filename, save_mkdir
from Tools.LinuxTools import CGAS


class Exonerate(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "exonerate", path=path, max_threads=max_threads)

    @staticmethod
    def parse_common_options(model, query_type=None, target_type=None,
                             show_alignment=None, show_sugar=True, show_cigar=None,
                             show_vulgar=None, show_query_gff=None, show_target_gff=None,
                             number_of_results_to_report=None, other_options=None):

        options = " --model %s" % model
        options += " --showalignment" if show_alignment else ""
        options += " --showsugar" if show_sugar else ""
        options += " --showcigar" if show_cigar else ""
        options += " --showvulgar" if show_vulgar else ""
        options += " --showquerygff" if show_query_gff else ""
        options += " --showtargetgff" if show_target_gff else ""
        options += " -Q %s" % query_type if query_type else ""
        options += " -T %s" % target_type if target_type else ""
        options += " -n %i" % number_of_results_to_report if number_of_results_to_report else ""
        options += " %s" % other_options if other_options else ""

        return options

    def parallel_alignment(self, query_file, target_file, model, num_of_recs_per_file=1000,
                           show_alignment=None, show_sugar=True, show_cigar=None,
                           show_vulgar=None, show_query_gff=None, show_target_gff=None,
                           store_intermediate_files=False,
                           splited_fasta_dir="splited_fasta_dir", splited_result_dir="splited_output",
                           number_of_results_to_report=None,
                           other_options=None,
                           converted_output_dir="converted_output"):
        splited_filename = split_filename(query_file)
        self.split_fasta(query_file, splited_fasta_dir, num_of_recs_per_file=num_of_recs_per_file,
                         output_prefix=splited_filename[1])

        common_options = self.parse_common_options(model, show_alignment=show_alignment,
                                                   show_sugar=show_sugar, show_cigar=show_cigar,
                                                   show_vulgar=show_vulgar, show_query_gff=show_query_gff,
                                                   show_target_gff=show_target_gff,
                                                   number_of_results_to_report=number_of_results_to_report,
                                                   other_options=other_options)

        options_list = []
        splited_files = os.listdir(splited_fasta_dir)

        save_mkdir(splited_result_dir)
        #save_mkdir(converted_output_dir)

        for filename in splited_files:
            filename_list = split_filename(filename)
            options = common_options
            options += " -q %s/%s" % (splited_fasta_dir, filename)
            options += " -t %s" % target_file
            options += " > %s/%s.output" % (splited_result_dir, filename_list[1])
            options_list.append(options)

        self.parallel_execute(options_list)

        """
        for filename in splited_files:

            trf_output_file = "%s/%s.%i.%i.%i.%i.%i.%i.%i.dat" % (splited_result_dir, filename,
                                                                  matching_weight, mismatching_penalty,
                                                                  indel_penalty, match_probability,
                                                                  indel_probability,
                                                                  min_alignment_score, max_period)

            self.convert_trf_report(trf_output_file, "%s/%s" % (converted_output_dir, filename))

        for suffix in (".rep", ".gff", ".simple.gff", ".short.tab", ".wide.tab"):
            file_str = ""
            merged_file = "%s%s" % (output_prefix, suffix)
            for filename in splited_files:
                file_str += " %s/%s%s" % (converted_output_dir, filename, suffix)
            CGAS.cat(file_str, merged_file)
        """
        if not store_intermediate_files:
            shutil.rmtree(splited_fasta_dir)
            #shutil.rmtree(splited_result_dir)
            #shutil.rmtree(converted_output_dir)

if __name__ == "__main__":
    pass
