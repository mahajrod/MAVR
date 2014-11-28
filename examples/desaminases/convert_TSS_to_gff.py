#!/usr/bin/env python

import os


tss_file = "/home/mahajrod/genetics/desaminases/data/TSS/TSS_sorted.t"
out_gff = "/home/mahajrod/genetics/desaminases/data/TSS/TSS_sorted.gff"
# sort -g -k 2 -k 4 TSS_no_header.t > TSS_no_header_sorted.t
source = "Zhenyu"
chr_names_dict = {
                  "1": "chrI",
                  "2": "chrII",
                  "3": "chrIII",
                  "4": "chrIV",
                  "5": "chrV",
                  "6": "chrVI",
                  "7": "chrVII",
                  "8": "chrVIII",
                  "9": "chrIX",
                  "10": "chrX",
                  "11": "chrXI",
                  "12": "chrXII",
                  "13": "chrXIII",
                  "14": "chrXIV",
                  "15": "chrXV",
                  "16": "chrXVI"
                  }

with open(tss_file, "r") as in_fd:
    with open(out_gff, "w") as out_fd:
        for line in in_fd:
            line_list = line.strip().split("\t")
            if len(line_list) == 1:
                break
            if line_list[0] == "ID":
                out_fd.write("##gff-version 3\n")
                continue

            TSS_name = [x.strip() for x in line_list[6].split(",")]
            id = TSS_name[0]

            TSS_common_name = [x.strip() for x in line_list[7].split(",")]
            out_fd.write("%s\t%s\t%s\t%s\t%s\t.\t%s\t.\tID=%s;Name=%s;Alias=%s;endConfidence=%s;Checked=%s\n" % (chr_names_dict[line_list[1]],
                                       source, line_list[5], line_list[3], line_list[4], line_list[2],
                                       id, ",".join(TSS_name), ",".join(TSS_common_name), line_list[8], line_list[9]
                                       ))

