#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool

from Parsers.TRF import CollectionTRF


class TRF(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "trf", path=path, max_threads=max_threads)

    def search_tandem_repeats(self, query_file, matching_weight=2, mismatching_penalty=7, indel_penalty=7,
                              match_probability=80, indel_probability=10, min_alignment_score=50, max_period=500,
                              report_flanking_sequences=False, make_dat_file=True, disable_html_output=True):

        options = " %s" % query_file
        options += " %i" % matching_weight
        options += " %i" % mismatching_penalty
        options += " %i" % indel_penalty
        options += " %i" % match_probability
        options += " %i" % indel_probability
        options += " %i" % min_alignment_score
        options += " %i" % max_period
        options += " -f" if report_flanking_sequences else ""
        options += " -d" if make_dat_file else ""
        options += " -h" if disable_html_output else ""

        self.execute(options)

    @staticmethod
    def convert_trf_report(trf_report, output_prefix):

        trf_collection = CollectionTRF(trf_file=trf_report, from_file=True)

        trf_collection.write("%s.rep" % output_prefix)
        trf_collection.write_gff("%s.gff" % output_prefix)
        trf_collection.write_short_table("%s.short.tab" % output_prefix)
        trf_collection.write_wide_table("%s.wide.tab" % output_prefix)

if __name__ == "__main__":
    pass
