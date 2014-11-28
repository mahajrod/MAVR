#!/usr/bin/env python
import numpy as np
from collections import OrderedDict


class TwoLvlDict(OrderedDict):

    def table_form(self, absent_symbol="0", sort=True):
        first_level_keys = list(self.keys())
        second_level_keys = set([])
        for fl_key in first_level_keys:
            for sl_key in self[fl_key]:
                second_level_keys.add(sl_key)
        if sort:
            first_level_keys.sort()
            second_level_keys = sorted(list(second_level_keys))
        string = ".\t" + "\t".join([str(x) for x in first_level_keys]) + "\n"

        for sl_key in second_level_keys:
            key_counts_list = []
            for fl_key in first_level_keys:
                if sl_key not in self[fl_key]:
                    key_counts_list.append(absent_symbol)
                else:
                    key_counts_list.append(str(self[fl_key][sl_key]))
            string += str(sl_key) + "\t" + "\t".join(key_counts_list) + "\n"
        return string

    def write(self, out_filename, absent_symbol="0"):
        with open(out_filename, "w") as out_fd:
            out_fd.write(self.table_form(absent_symbol=absent_symbol))

