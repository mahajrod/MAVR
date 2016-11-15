#!/usr/bin/env python
__author__ = 'mahajrod'
from collections import OrderedDict


class CoockiecutterReport():
    def __init__(self, coockiecutter_report_file):
        self.stat = OrderedDict()
        self.stat["left"] = OrderedDict()
        self.stat["right"] = OrderedDict()

        with open(coockiecutter_report_file, "r") as in_fd:
            in_fd.readline()
            in_fd.readline()
            for i in range(0, 6):
                tmp = in_fd.readline().strip().split()
                if i == 3:
                    self.stat["left"][tmp[0]] = float(tmp[1])
                else:
                    self.stat["left"][tmp[0]] = int(tmp[1])
            in_fd.readline()
            for i in range(0, 6):
                tmp = in_fd.readline().strip().split()
                if i == 3:
                    self.stat["right"][tmp[0]] = float(tmp[1])
                else:
                    self.stat["right"][tmp[0]] = int(tmp[1])

        self.input_pairs = self.stat["left"]["ok"] + self.stat["left"]["adapter"] + self.stat["left"]["n"]
        self.retained_pairs = self.stat["left"]["pe"]
