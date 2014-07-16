#!/usr/bin/env python
__author__ = 'mahajrod'

#TODO: write parser


class FastQCRecord():
    def __init__(self, fastqc_report_file):
        with open(fastqc_report_file, "r") as in_fd:
            self.version = in_fd.readline().strip().split("\t")[1]
            self.report = {}
            for line in in_fd:
                temp = line.strip()
                if temp[0:2] == ">>":
                    module_name, quality = temp[2:].split("\t")
                    self.report[module_name] = {"quality": quality, "header": [], "data": []}
                    for line in in_fd:
                        temp = line.strip()
                        if temp == ">>END_MODULE":
                            break
                        if temp[0] == "#":
                            self.report[module_name]["header"].append(temp[1:].split("\t"))
                        else:
                            self.report[module_name]["data"].append(temp.split("\t"))

if __name__ == "__main__":
    FastQC_file = FastQCRecord("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all/N080-LAN210-Can-A1-NA-RUN6-D6/N080-LAN210-Can-A1-NA-RUN6-D6-_TCCTGAGC-TATCCTCT_L002_R2_001_fastqc/fastqc_data.txt")
    print(FastQC_file.report["Basic Statistics"])