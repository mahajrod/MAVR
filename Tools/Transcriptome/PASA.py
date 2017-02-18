#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from Tools.Abstract import Tool


class PASA(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "Launch_PASA_pipeline.pl", path=path, max_threads=max_threads)

    def seqclean(self, fasta_with_transcripts, output_file=None, report_file=None, seq_number_per_search_slice=None,
                 seq_to_trim_file_list=None, disable_trimming_N_rich_ends=None,
                 disable_trashing_of_low_quality_sequences=None,  disable_trimming_polyAT_tails=None,
                 reject_sequences_shorter=None, contamination_sequences_file_list=None):
        """
        """
        options = " %s" % fasta_with_transcripts
        options += " -n %i" % seq_number_per_search_slice if seq_number_per_search_slice else ""
        options += (" -v %s" % ",".join(seq_to_trim_file_list)) if seq_to_trim_file_list else ""
        options += " -l %i" % reject_sequences_shorter if reject_sequences_shorter else ""
        options += (" -s %s" ",".join(contamination_sequences_file_list)) if contamination_sequences_file_list else ""

        options += " -r %s" % report_file if report_file else ""
        options += " -o %s" % output_file if output_file else ""

        options += " -N" if disable_trimming_N_rich_ends else ""
        options += " -M" if disable_trashing_of_low_quality_sequences else ""
        options += " -A" if disable_trimming_polyAT_tails else ""

        self.execute(options, cmd="seqclean")

if __name__ == "__main__":
    pass
