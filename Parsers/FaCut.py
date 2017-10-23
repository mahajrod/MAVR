#!/usr/bin/env python
__author__ = 'mahajrod'
import numpy as np
from collections import OrderedDict

from Routines import MatplotlibRoutines

from CustomCollections.GeneralCollections import TwoLvlDict


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
                                      "both_discarded": 8,
                                      "total": 9})
        self.tile_table = []
        self.machine_id_list = []
        self.flowcell_id_list = []

        with open(facut_report_file, "r") as in_fd:
            tmp = ["", ""]
            while tmp[0] != "Paires retained:":
                tmp = in_fd.readline().strip().split("\t")

            self.retained_pairs = int(tmp[1])

            self.retained_forward_only = int(in_fd.readline().strip().split("\t")[1])
            self.retained_reverse_only = int(in_fd.readline().strip().split("\t")[1])
            self.both_discarded = int(in_fd.readline().strip().split("\t")[1])
            self.input_pairs = self.retained_pairs + self.retained_forward_only + self.retained_reverse_only + self.both_discarded
            self.retained_pairs_fraction = float(self.retained_pairs) / float(self.input_pairs)

            tmp = "aaaaaaaaaaaaaaaaabbbbbbbbbbbbb"
            while tmp[:13] != "instrument_id":
                tmp = in_fd.readline()

            for line in in_fd:
                tmp = line.strip().split("\t")
                for i in (1, 3, 4, 5, 6, 7, 8):
                    tmp[i] = int(tmp[i])
                tmp.append(tmp[5] + tmp[6] + tmp[7] + tmp[8])
                if tmp[0] not in self.machine_id_list:
                    self.machine_id_list.append(tmp[0])
                tmp[0] = self.machine_id_list.index(tmp[0])  # replace machine id by index

                if tmp[2] not in self.flowcell_id_list:
                    self.flowcell_id_list.append(tmp[2])
                tmp[2] = self.flowcell_id_list.index(tmp[2])  # replace flowcell id by index

                self.tile_table.append(tmp)
        #print self.tile_table
        self.tile_table = np.array(self.tile_table)

        self.minimum_retained_pairs_in_tiles_fraction = min(self.get_fraction_of_retained_pairs_per_tile())

        self.lane_table = self.get_lane_table()
        self.full_lane_id_list = self.get_full_lane_id_list()
        self.short_lane_id_list = self.get_short_lane_id_list() # WARNING! only lane number is stored here, so if there are lanes with same numbers of several runs they will collaps

        self.more_than_one_machine = True if len(self.machine_id_list) > 1 else False
        self.more_than_one_flowcell = True if len(self.flowcell_id_list) > 1 else False
        self.more_than_one_lane = True if len(self.lane_table) > 1 else False

    def get_fraction_of_retained_pairs_per_tile(self):
        return np.array(self.tile_table[:, 5], dtype=float) / np.array(self.tile_table[:, 9], dtype=float)

    def draw_fraction_of_retained_pairs_per_tile_histogram(self, output_prefix):
        data = self.get_fraction_of_retained_pairs_per_tile()
        MatplotlibRoutines.percent_histogram(data, output_prefix,
                                             n_bins=20, title="Distribution of retained pairs per tile",
                                             xlabel="Fraction of retained pairs", ylabel="Number of tiles", label=None,
                                             extensions=("png", "svg"), legend=None, legend_location="best",
                                             input_mode="fraction", xmax=None, xmin=None)

    def get_lane_table(self):
        # split tile_table by lanes
        lane_tile_list = np.split(self.tile_table,
                                  np.unique(np.where(np.diff(self.tile_table[:, 0:self.table_ids["tile"]], axis=0))[0]) + 1)
        """
        print len((self.tile_table[:, 0:self.table_ids["tile"]]))
        print list(np.where(np.diff(self.tile_table[:, 0:self.table_ids["tile"]], axis=0)))
        print "CCCCCCCCCCCCCCCCCCCCCC"

        print list(np.diff(self.tile_table[:, 0:self.table_ids["tile"]], axis=0))
        print "BBBBBBBBBBBBBBB"


        print lane_tile_list
        print "AAAAAAAAAAAAAAAAAAAAAAAAAAa"
        print "\n"
        """
        #print lane_tile_list
        lane_table = []
        for lane_tile_table in lane_tile_list:
            #print lane_tile_table
            #print "\n"
            lane_ids = list(lane_tile_table[0][:self.table_ids["tile"]])
            lane_values = list(np.sum(lane_tile_table[:, self.table_ids["both_retained"]:], axis=0))

            lane_table.append(lane_ids + [0] + lane_values)

        return np.array(lane_table)

    def lane_table_string(self):
        lane_str = "\t".join(self.table_ids.keys()) + "\n"
        for lane_line in self.lane_table:
            tmp = list(lane_line)
            tmp[0] = self.machine_id_list[tmp[0]]
            tmp[2] = self.flowcell_id_list[tmp[2]]
            lane_str += "\t".join(map(str, tmp)) + "\n"

        return lane_str

    def get_full_lane_id_list(self):
        lane_id_list = []
        for lane_line in self.lane_table:
            tmp = list(lane_line)
            lane_id_list.append(":".join([self.machine_id_list[tmp[0]],
                                          str(tmp[1]),
                                          self.flowcell_id_list[tmp[2]],
                                          str(tmp[3])]))
        return lane_id_list

    def get_short_lane_id_list(self):
        return list(np.unique(self.tile_table[:, self.table_ids["lane_number"]]))


class FaCutReportCollection(OrderedDict):
    def __init__(self, report_dict=None, file_dict=None):
        OrderedDict.__init__(self)
        if report_dict:
            for report_id in report_dict:
                self[report_id] = report_dict
        if file_dict:
            for file_id in file_dict:
                if file_id in self:
                    raise KeyError("Report with same id(%s) is already present in FaCutReportCollection" % str(file_id))
                self[file_id] = FaCutReport(file_dict[file_id])

    def get_general_stats(self):
        stat_dict = TwoLvlDict()

        for report_id in self:
            stat_dict[report_id] = OrderedDict()
            stat_dict[report_id]["machine_number"] = len(self[report_id].machine_id_list)
            stat_dict[report_id]["machine_ids"] = self[report_id].machine_id_list
            stat_dict[report_id]["flowcell_number"] = len(self[report_id].flowcell_id_list)
            stat_dict[report_id]["flowcell_ids"] = self[report_id].flowcell_id_list
            stat_dict[report_id]["lane_number"] = len(self[report_id].lane_table)
            stat_dict[report_id]["full_lane_ids"] = self[report_id].full_lane_id_list
            stat_dict[report_id]["short_lane_ids"] = self[report_id].short_lane_id_list
            stat_dict[report_id]["input_pairs"] = self[report_id].input_pairs
            stat_dict[report_id]["retained_pairs"] = self[report_id].retained_pairs
            stat_dict[report_id]["retained_pairs_fraction"] = self[report_id].retained_pairs_fraction
            stat_dict[report_id]["retained_forward_only"] = self[report_id].retained_forward_only
            stat_dict[report_id]["retained_reverse_only"] = self[report_id].retained_reverse_only
            stat_dict[report_id]["both_discarded"] = self[report_id].both_discarded
            stat_dict[report_id]["min_retained_pairs_in_tiles_fraction"] = self[report_id].minimum_retained_pairs_in_tiles_fraction

        return stat_dict

    def write_general_stats(self, output_file):
        stat_dict = self.get_general_stats()
        stat_dict.write(output_file, absent_symbol=".", sort=False)
