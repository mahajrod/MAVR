#!/usr/bin/env python
__author__ = 'mahajrod'
from collections import OrderedDict


#TODO: think about parsing data in modules - at moment all data are counted as strings
#tested with FastQC v0.10.1
class FastQCModule():
    def __init__(self, name, status, header_list, data, info_list=None):
        self.name = name
        self.status = status
        self.header_list = header_list
        self.data = data
        self.info_list = info_list

    def __str__(self):
        string = ">>%s\t%s\n" % (self.name, self.status)
        string += ("#" + "\t".join(self.info_list) + "\n") if self.info_list else ""
        string += "#" + "\t".join(self.header_list) + "\n"
        keys = list(self.data.keys())
        for i in range(0, len(self.data[keys[0]])):
            entry = []
            for key in keys:
                entry.append(str(self.data[key][i]))
            string += "\t".join(entry) + "\n"

        string += ">>END_MODULE\n"
        return string


class FastQCReport():
    def __init__(self, fastqc_report_file):
        with open(fastqc_report_file, "r") as in_fd:
            self.version = in_fd.readline().strip().split("\t")[1]
            self.report = OrderedDict({})
            for line in in_fd:
                temp = line.strip()
                if temp[0:2] == ">>":
                    module_name, status = temp[2:].split("\t")
                    self.report[module_name] = self.init_module(in_fd, module_name, status)

    def init_module(self, fd, name, status):
        info_list = fd.readline().strip()[1:].split("\t") if name == "Sequence Duplication Levels" else None
        header_list = fd.readline().strip()[1:].split("\t")
        data = OrderedDict({})
        for entry in header_list:
            data[entry] = []
        for line in fd:
            if line[:2] == ">>":
                return FastQCModule(name, status, header_list, data, info_list=info_list)
            temp = line.strip().split("\t")
            for (key, value) in zip(header_list, temp):
                data[key].append(value)

    def overrepresented_sequences(self):
        return self.report["Overrepresented sequences"].data["Sequence"]

    def overrepresented_kmers(self):
        return self.report["Kmer Content"].data["Sequence"]

if __name__ == "__main__":
    FastQC_file = FastQCReport("../example_data/FastQC/fastqc_data.txt")
    print(FastQC_file.overrepresented_sequences())
    print(FastQC_file.overrepresented_kmers())