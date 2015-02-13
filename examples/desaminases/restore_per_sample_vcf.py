import os

from Tools.Bedtools import  Intersect
from  Parsers.VCF import CollectionVCF, HeaderVCF, RecordVCF



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

HAP_list = ["N006-LAN211-Can-HAP-NA-RUN1",
            "N007-LAN211-Can-HAP-NA-RUN1",
            "N014-LAN211-Can-HAP-NA-RUN2",
            "N015-LAN211-Can-HAP-NA-RUN2",
            "N016-LAN211-Can-HAP-NA-RUN2",
            "N017-LAN211-Can-HAP-NA-RUN2",
            "N018-LAN211-Can-HAP-NA-RUN2",
            "N019-LAN211-Can-HAP-NA-RUN2",
            "N050-LAN211-Can-HAP-NA-RUN4",
            "N051-LAN211-Can-HAP-NA-RUN4"
            ]

HAP_sub1_list = ["N081-LAN210_sub1KanMX-Can-HAP-NA-RUN6",
                 "N082-LAN210_sub1KanMX-Can-HAP-NA-RUN6"
                 ]

sample_set_dict = {"PmCDA1_3d": PmCDA1_3d_list,
                   "PmCDA1_6d": PmCDA1_6d_list,
                   "A1_3d": A1_3d_list,
                   "A1_6d": A1_6d_list,
                   "A3G_3d": A3G_3d_list,
                   "AID_3d": AID_3d_list,
                   "AID_6d": AID_6d_list,
                   "PmCDA1_sub1_3d": PmCDA1_sub1_3d_list,
                   "PmCDA1_sub1_6d": PmCDA1_sub1_6d_list,
                   "HAP": HAP_list,
                   "HAP_sub1": HAP_sub1_list
                   }

sample_set_names_list = ["PmCDA1_3d",
                         "PmCDA1_6d",
                         "A1_3d",
                         "A1_6d",
                         "A3G_3d",
                         "AID_3d",
                         "AID_6d",
                         "PmCDA1_sub1_3d",
                         "PmCDA1_sub1_6d",
                         "HAP",
                         "HAP_sub1"
                         ]


combined_vcf_suffix = "_adjusted_cluster_mutations.vcf"
combined_3_vcf_suffix = "_adjusted_3+_cluster_mutations.vcf"

samples_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all/"
samples_subdir = "alignment_LAN210_v0.10m/"
sample_suffix = "_GATK_best_SNP.vcf"

workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"
combined_subdir = "clustering/"


for sample_set in sample_set_names_list:
    os.chdir(workdir)
    os.chdir(sample_set)
    os.system("mkdir -p per_sample_vcf")
    sample_set_mutations = CollectionVCF(from_file=True, vcf_file=combined_subdir + sample_set + combined_vcf_suffix)
    #sample_set_3_mutations = CollectionVCF(from_file=True, vcf_file=combined_subdir + sample_set + combined_3_vcf_suffix)
    per_sample_mutations = {}
    samples_list = sample_set_mutations.samples
    print("Handling %s" % sample_set)
    for sample in sample_set_mutations.samples:
        print("Handling %s" % sample)
        sample_mutations = CollectionVCF(from_file=True, vcf_file=samples_dir + sample + "/" + samples_subdir  + sample + sample_suffix)
        per_sample_mutations[sample] = CollectionVCF(metadata=sample_mutations.metadata, record_list=None,
                                                     header=sample_mutations.header,
                                                     vcf_file=None, samples=None, from_file=False, external_metadata=None)

        sample_index = samples_list.index(sample_mutations.samples[0])
        for sample_record in sample_mutations:
            for sample_set_record in sample_set_mutations:
                if sample_record.chrom != sample_set_record.chrom or sample_set_record.pos != sample_record.pos:
                    continue
                if sample_set_record.samples_list[sample_index]["GT"] != "./.":
                    per_sample_mutations[sample].records.append(sample_record)

        per_sample_mutations[sample].write(workdir + sample_set  + "/per_sample_vcf/" + sample + ".vcf")


