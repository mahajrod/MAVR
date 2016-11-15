#!/usr/bin/env python
__author__ = 'mahajrod'
import numpy as np
from collections import OrderedDict

from Routines import MatplotlibRoutines


class FaCutReport:
    def __init__(self, facut_report_file):
        self.table_ids = OrderedDict({"instrument_id": 0,
                                      "run_number": 1,
                                      "flowcell_id": 2,
                                      "lane_number": 3,
                                      "tile": 4,
                                      "both_retained": 5,
                                      "forward_only": 6,
                                      "reverse_only": 7,
                                      "both_discarded": 8})
        self.tile_table = []
        with open(facut_report_file, "r") as in_fd:
            tmp = ["", ""]
            while tmp[0] != "Paires retained:":
                tmp = in_fd.readline().strip().split("\t")

            self.retained_pairs = int(tmp[1])

            self.retained_forward_only = int(in_fd.readline().strip().split("\t")[1])
            self.retained_reverse_only = int(in_fd.readline().strip().split("\t")[1])
            self.both_discarded = int(in_fd.readline().strip().split("\t")[1])
            self.input_pairs = self.retained_pairs + self.retained_forward_only + self.retained_reverse_only + self.both_discarded
            tmp = "aaaaaaaaaaaaaaaaabbbbbbbbbbbbb"
            while tmp[:13] != "instrument_id":
                tmp = in_fd.readline()

            for line in in_fd:
                tmp = line.strip().split("\t")
                for i in (1, 3, 5, 6, 7, 8):
                    tmp[i] = int(tmp[i])
                tmp.append(tmp[5] + tmp[6] + tmp[7] + tmp[8])

                self.tile_table.append(tuple(tmp))

        self.tile_table = np.array(self.tile_table, dtype='|S50, u4, S50, u4, S6, u8, u8, u8, u8, u8')
        print self.tile_table
        print self.tile_table[1]
        print self.tile_table[1][5]

        print self.tile_table[:][5]

    def get_fraction_of_retained_pairs_per_tile(self):
        return self.tile_table[:, 5] / self.tile_table[:, 9]

    def draw_fraction_of_retained_pairs_per_tile_histogram(self, output_prefix):
        MatplotlibRoutines.percent_histogram(self.get_fraction_of_retained_pairs_per_tile(), output_prefix,
                                             n_bins=20, title="Distribution of retained pairs per pile",
                                             xlabel="Fraction", ylabel="Number", label=None,
                                             extensions=("png", "svg"), legend=None, legend_location="best",
                                             input_mode="fraction", xmax=None, xmin=None)

