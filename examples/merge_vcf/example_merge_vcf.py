#!/usr/bin/env python
import os
from Tools.GATK import CombineVariants


PmCDA1_3d_list = ["N010-LAN210-Can-PmCDA1-NA-RUN2-D3",
                  "N011-LAN210-Can-PmCDA1-NA-RUN2-D3",
                  "N012-LAN210-Can-PmCDA1-NA-RUN2-D3",
                  "N013-LAN210-Can-PmCDA1-NA-RUN2-D3",
                  "N061-LAN210-Can-PmCDA1-NA-RUN5-D3",
                  "N062-LAN210-Can-PmCDA1-NA-RUN5-D3",
                  "N065-LAN210-Can-PmCDA1-NA-RUN6-D3",
                  "N066-LAN210-Can-PmCDA1-NA-RUN6-D3",
                  "N058-LAN210-Can-PmCDA1-Feb13-RUN5-D3",
                  "N059-LAN210-FOA-PmCDA1-Feb13-RUN5-D3",
                  "N060-LAN210-FOA-PmCDA1-Feb13-RUN5-D3"
                  ]

PmCDA1_6d_list = ["N070-LAN210-Can-PmCDA1-NA-RUN6-D6",
                  "N071-LAN210-Can-PmCDA1-NA-RUN6-D6",
                  "N072-LAN210-Can-PmCDA1-NA-RUN6-D6",
                  "N073-LAN210-Can-PmCDA1-NA-RUN6-D6"
                  ]

A1_3d_list = ["N034-LAN210-Can-A1-Oct12-RUN4-D3",
              "N035-LAN210-Can-A1-Oct12-RUN4-D3",
              "N036-LAN210-FOA-A1-Oct12-RUN4-D3",
              "N037-LAN210-FOA-A1-Oct12-RUN4-D3"
              ]

A1_6d_list = ["N079-LAN210-Can-A1-NA-RUN6-D6",
              "N080-LAN210-Can-A1-NA-RUN6-D6"
              ]


A3G_3d_list = ["N038-LAN210-Can-A3G-Oct12-RUN4-D3",
               "N039-LAN210-Can-A3G-Oct12-RUN4-D3"
               ]

AID_3d_list = ["N040-LAN210-Can-AID-Oct12-RUN4-D3",
               "N041-LAN210-Can-AID-Oct12-RUN4-D3"
               ]

AID_6d_list = ["N077-LAN210-Can-AID-NA-RUN6-D6",
               "N078-LAN210-Can-AID-NA-RUN6-D6"
               ]

PmCDA1_sub1_3d_list = ["N067-LAN210_sub1KanMX-Can-PmCDA1-NA-RUN6-D3",
                       "N068-LAN210_sub1KanMX-Can-PmCDA1-NA-RUN6-D3",
                       "N069-LAN210_sub1KanMX-Can-PmCDA1-NA-RUN6-D3",
                        ]

PmCDA1_sub1_6d_list = ["N074-LAN210_sub1KanMX-Can-PmCDA1-NA-RUN6-D6",
                       "N075-LAN210_sub1KanMX-Can-PmCDA1-NA-RUN6-D6",
                       "N076-LAN210_sub1KanMX-Can-PmCDA1-NA-RUN6-D6"
                       ]

sample_set_list = [PmCDA1_3d_list,
                   PmCDA1_6d_list,
                   A1_3d_list,
                   A1_6d_list,
                   A3G_3d_list,
                   AID_3d_list,
                   AID_6d_list,
                   PmCDA1_sub1_3d_list,
                   PmCDA1_sub1_6d_list]
sample_set_names_list = ["PmCDA1_3d",
                         "PmCDA1_6d",
                         "A1_3d",
                         "A1_6d",
                         "A3G_3d",
                         "AID_3d",
                         "AID_6d",
                         "PmCDA1_sub1_3d",
                         "PmCDA1_sub1_6d"]

samples_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all/"
subdir = "alignment_LAN210_v0.10m/"
suffix = "_GATK_best_merged.vcf"
result_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"
gatk_dir = "/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0/"
reference = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta"
os.chdir(result_dir)

for samples_set, sample_set_name in zip(sample_set_list, sample_set_names_list):
    vcf_list = ["%s%s/%s%s%s" % (samples_dir, sample, subdir, sample, suffix) for sample in samples_set]
    print("Combining:\n\t" + "\n\t".join(vcf_list))
    CombineVariants.combine_from_same_source(gatk_dir, reference, vcf_list, sample_set_name + ".vcf", num_threads=4)

