#!/usr/bin/env python
__author__ = 'mahajrod'
from collections import OrderedDict


# adapter for 2 versions of Coockiecutter
class CoockiecutterReport():
    def __init__(self, coockiecutter_report_file):
        self.stat = OrderedDict()
        self.stat["left"] = OrderedDict()
        self.stat["right"] = OrderedDict()

        with open(coockiecutter_report_file, "r") as in_fd:
            tmp = ""

            for read_position in "left", "right":
                while "ok" not in tmp:
                    tmp = in_fd.readline()

                while "\t" in tmp:
                    tmp = tmp.strip().split("\t")
                    if "%" in tmp[1]:
                        self.stat[read_position][tmp[0]] = float(tmp[1].replace("%", ""))
                    else:
                        self.stat[read_position][tmp[0]] = int(tmp[1])
                    tmp = in_fd.readline()

        self.match_key = "match" if "match" in self.stat["left"] else "adapter"
        self.retained_pairs_key = "pe" if "pe" in self.stat["left"] else "paired-end reads"
        self.input_pairs = self.stat["left"]["ok"] + self.stat["left"][self.match_key] + self.stat["left"]["n"]
        self.retained_pairs = self.stat["left"][self.retained_pairs_key]
