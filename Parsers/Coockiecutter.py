#!/usr/bin/env python
__author__ = 'mahajrod'
from collections import OrderedDict


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

        self.input_pairs = self.stat["left"]["ok"] + self.stat["left"]["adapter"] + self.stat["left"]["n"]
        self.retained_pairs = self.stat["left"]["pe"]
