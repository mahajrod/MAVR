#!/usr/bin/env python
__author__ = 'mahajrod'
from collections import OrderedDict
from CustomCollections.GeneralCollections import TwoLvlDict


class KrATERReport(OrderedDict):
    def __init__(self, krater_report_file):
        OrderedDict.__init__(self)
        #self = OrderedDict()

        with open(krater_report_file, "r") as in_fd:
            for line in in_fd:
                tmp = line.strip().split("\t")
                self[tmp[0]] = tmp[1]

        self["Number of distinct kmers"] = int(self["Number of distinct kmers"])
        self["Number of distinct kmers with errors"] = int(self["Number of distinct kmers with errors"])
        self["Fraction of distinct kmers with errors"] = float(self["Fraction of distinct kmers with errors"])

        self["Total number of kmers"] = int(self["Total number of kmers"])
        self["Total number of kmers with errors"] = int(self["Total number of kmers with errors"])
        self["Fraction of kmers with errors"] = float(self["Fraction of kmers with errors"])

        self["Width of first peak"] = int(self["Width of first peak"])
        self["Mean kmer multiplicity in first peak"] = float(self["Mean kmer multiplicity in first peak"])
        self["Kmer multiplicity at first maximum"] = float(self["Kmer multiplicity at first maximum"])
        self["Standard deviation of kmer multiplicity in first peak"] = float(self["Standard deviation of kmer multiplicity in first peak"])
        self["Variance coefficient of kmer multiplicity in first peak"] = float(self["Variance coefficient of kmer multiplicity in first peak"])
        if "Estimated genome size, bp" in self: # for compatibility with older versions that don't count genome size
            self["Estimated genome size, bp"] = float(self["Estimated genome size, bp"]) if self["Estimated genome size, bp"] != 'NA' else None


class KrATERReportCollection(OrderedDict):
    def __init__(self, report_dict=None, file_dict=None):
        OrderedDict.__init__(self)
        if report_dict:
            for report_id in report_dict:
                self[report_id] = report_dict
        if file_dict:
            for file_id in file_dict:
                if file_id in self:
                    raise KeyError("Report with same id(%s) is already present in KrATERReportCollection" % str(file_id))
                self[file_id] = KrATERReport(file_dict[file_id])

    def get_general_stats(self):
        stat_dict = TwoLvlDict()

        for report_id in self:
            stat_dict[report_id] = OrderedDict()

            stat_dict[report_id]["Number of distinct kmers"] = self[report_id]["Number of distinct kmers"]
            stat_dict[report_id]["Number of distinct kmers"] = self[report_id]["Number of distinct kmers"]
            stat_dict[report_id]["Fraction of distinct kmers with errors"] = self[report_id]["Fraction of distinct kmers with errors"]

            stat_dict[report_id]["Total number of kmers"] = self[report_id]["Total number of kmers"]
            stat_dict[report_id]["Total number of kmers with errors"] = self[report_id]["Total number of kmers with errors"]
            stat_dict[report_id]["Fraction of kmers with errors"] = self[report_id]["Fraction of kmers with errors"]

            stat_dict[report_id]["Width of first peak"] = self[report_id]["Width of first peak"]
            stat_dict[report_id]["Mean kmer multiplicity in first peak"] = self[report_id]["Mean kmer multiplicity in first peak"]
            stat_dict[report_id]["Kmer multiplicity at first maximum "] = self[report_id]["Kmer multiplicity at first maximum"]
            stat_dict[report_id]["Standard deviation of kmer multiplicity in first peak"] = self[report_id]["Standard deviation of kmer multiplicity in first peak"]
            stat_dict[report_id]["Variance coefficient of kmer multiplicity in first peak"] = self[report_id]["Variance coefficient of kmer multiplicity in first peak"]
            if "Estimated genome size, bp" in self[report_id]:
                stat_dict[report_id]["Estimated genome size,bp"] = self[report_id]["Estimated genome size, bp"]

        return stat_dict

    def write_general_stats(self, output_file):
        stat_dict = self.get_general_stats()
        stat_dict.write(output_file, absent_symbol=".", sort=False)
