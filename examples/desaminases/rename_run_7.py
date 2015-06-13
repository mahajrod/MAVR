#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

import os

from collections import  OrderedDict
samples_dict = {"Sample_1": "N085-LAN210-Can-PmCDA1-NA-RUN7-D1",
                "Sample_2": "N086-LAN210-Can-PmCDA1-NA-RUN7-D1",
                "Sample_3": "N087-LAN210-Can-PmCDA1-NA-RUN7-D1",
                "Sample_4": "N088-LAN210-Can-PmCDA1-NA-RUN7-D1",
                "Sample_5": "N089-LAN210-Can-PmCDA1-NA-RUN7-D1",
                "Sample_6": "N090-LAN210-Can-PmCDA1-NA-RUN7-D1",
                "Sample_7": "N091-LAN210-Can-PmCDA1-NA-RUN7-D3",
                "Sample_8": "N092-LAN210-Can-PmCDA1-NA-RUN7-D6",
                "Sample_9": "N093-LAN210-Can-PmCDA1-NA-RUN7-D6",
                "Sample_10": "N094-LAN210-Can-PmCDA1-NA-RUN7-D6",
                "Sample_11": "N095-LAN210-Can-PmCDA1-NA-RUN7-D6",
                "Sample_12": "N096-LAN210-Can-PmCDA1-NA-RUN7-D6",
                "Sample_13": "N097-LAN210-Can-AID-NA-RUN7-D1",
                "Sample_14": "N098-LAN210-Can-AID-NA-RUN7-D1",
                #"Sample_15": "N099-LAN210-Can-AID-NA-RUN7-D1",
                "Sample_16": "N100-LAN210-Can-AID-NA-RUN7-D3",
                "Sample_17": "N101-LAN210-Can-AID-NA-RUN7-D3",
                "Sample_18": "N102-LAN210-Can-AID-NA-RUN7-D3",
                "Sample_19": "N103-LAN210-Can-AID-NA-RUN7-D6",
                "Sample_20": "N104-LAN210-Can-AID-NA-RUN7-D6",
                "Sample_21": "N105-LAN210-Can-AID-NA-RUN7-D6",
                "Sample_22": "N106-LAN210-Can-AID-NA-RUN7-D6",
                "Sample_23": "N107-LAN210-Can-A1-NA-RUN7-D1",
                "Sample_24": "N108-LAN210-Can-A1-NA-RUN7-D1",
                "Sample_25": "N109-LAN210-Can-A1-NA-RUN7-D1",
                "Sample_26": "N110-LAN210-Can-A1-NA-RUN7-D1",
                "Sample_27": "N111-LAN210-Can-A1-NA-RUN7-D1",
                "Sample_28": "N112-LAN210-Can-A1-NA-RUN7-D1",
                "Sample_29": "N113-LAN210-Can-A1-NA-RUN7-D3",
                "Sample_30": "N114-LAN210-Can-A1-NA-RUN7-D3",
                "Sample_31": "N115-LAN210-Can-A1-NA-RUN7-D3",
                "Sample_32": "N116-LAN210-Can-A1-NA-RUN7-D3",
                "Sample_33": "N117-LAN210-Can-A1-NA-RUN7-D3",
                "Sample_34": "N118-LAN210-Can-A1-NA-RUN7-D3",
                #"Sample_35": "N119-LAN210-Can-A1-NA-RUN7-D3",
                "Sample_36": "N120-LAN210-Can-A1-NA-RUN7-D6",
                "Sample_37": "N121-LAN210-Can-A1-NA-RUN7-D6",
                "Sample_38": "N122-LAN210-Can-A1-NA-RUN7-D6",
                "Sample_39": "N123-LAN210-Can-A1-NA-RUN7-D6",
                "Sample_40": "N124-LAN210-Can-A1-NA-RUN7-D6"}

raw_dir = "raw_names/"
workdir = "./"

for sample in samples_dict:
    print("Handling %s..." % sample)
    os.system("mkdir -p %s" % samples_dict[sample])
    sample_files = os.listdir(raw_dir + sample)
    for filename in sample_files:
        if "R1" in filename:
            os.system("cp %s%s/%s %s/%s_1.fastq" % (raw_dir, sample, filename, samples_dict[sample], samples_dict[sample]))
        elif "R2" in filename:
            os.system("cp %s%s/%s %s/%s_2.fastq" % (raw_dir, sample, filename, samples_dict[sample], samples_dict[sample]))
        else:
            os.system("cp %s%s/%s %s/" % (raw_dir, sample, filename, samples_dict[sample]))