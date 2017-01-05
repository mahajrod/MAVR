#!/usr/bin/env python
__author__ = 'mahajrod'
from collections import OrderedDict
from CustomCollections.GeneralCollections import TwoLvlDict

# adapter for 2 versions of Coockiecutter
class CoockiecutterReport(OrderedDict):
    def __init__(self, coockiecutter_report_file):
        #self = OrderedDict()
        OrderedDict.__init__(self)
        self["left"] = OrderedDict()
        self["right"] = OrderedDict()

        with open(coockiecutter_report_file, "r") as in_fd:
            tmp = ""

            for read_position in "left", "right":
                while "ok" not in tmp:
                    tmp = in_fd.readline()

                while "\t" in tmp:
                    tmp = tmp.strip().split("\t")
                    if "%" in tmp[1]:
                        self[read_position][tmp[0]] = float(tmp[1].replace("%", ""))
                    else:
                        self[read_position][tmp[0]] = int(tmp[1])
                    tmp = in_fd.readline()

        self.match_key = "match" if "match" in self["left"] else "adapter"
        self.retained_pairs_key = "pe" if "pe" in self["left"] else "paired-end reads"
        self.input_pairs = self["left"]["ok"] + self["left"][self.match_key] + self["left"]["n"]
        self.retained_pairs = self["left"][self.retained_pairs_key]
        self.retained_pairs_fraction = float(self.retained_pairs) / float(self.input_pairs)


class CoockiecutterReportCollection(OrderedDict):
    def __init__(self, report_dict=None, file_dict=None):
        OrderedDict.__init__(self)
        if report_dict:
            for report_id in report_dict:
                self[report_id] = report_dict
        if file_dict:
            for file_id in file_dict:
                if file_id in self:
                    raise KeyError("Report with same id(%s) is already present in CoockiecutterReportCollection" % str(file_id))
                self[file_id] = CoockiecutterReport(file_dict[file_id])

    def get_general_stats(self):
        stat_dict = TwoLvlDict()

        for report_id in self:
            stat_dict[report_id] = OrderedDict()

            stat_dict[report_id]["input_pairs"] = self[report_id].input_pairs
            stat_dict[report_id]["pairs_without_adapters"] = self[report_id].retained_pairs
            stat_dict[report_id]["pairs_without_adapters_fraction"] = self[report_id].retained_pairs_fraction

        return stat_dict

    def write_general_stats(self, output_file):
        stat_dict = self.get_general_stats()
        stat_dict.write(output_file, absent_symbol=".", sort=False)
