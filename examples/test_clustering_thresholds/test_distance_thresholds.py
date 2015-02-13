__author__ = 'mahajrod'
import os
from Parsers.VCF import CollectionVCF

if __name__ == "__main__":
    """
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all"

    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))

    suffix = "_GATK_best_merged.vcf"

    for sample in samples_list:
        print("Handling %s" % sample)

        os.chdir(workdir)
        os.chdir(sample)
        if "alignment_LAN210_v0.9m" not in os.listdir("."):
            continue
        os.chdir("alignment_LAN210_v0.9m")

        mutations = CollectionVCF(vcf_file=sample + suffix,
                                  from_file=True)
        print("Totaly %s mutations" % len(mutations))

        mutations.test_thresholds(extracting_method='distance', threshold=(50, 5000, 100),
                                  testing_dir="testing_threshold")
    """
    sample_set_names_list = ["PmCDA1_3d",
                             "HAP",
                             "PmCDA1_sub1_3d",
                             "PmCDA1_6d",
                             "HAP_sub1",
                             "PmCDA1_sub1_6d",
                             "A1_3d",
                             "A1_6d",
                             "A3G_3d",
                             "AID_3d",
                             "AID_6d"
                             ]
    suffix = "_SNP.vcf"
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/test_thresholds/"
    sample_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/SNP_annotated_raw_vcf/"

    for sample_set in sample_set_names_list:
        print("Handling %s" % sample_set)
        os.chdir(workdir)
        os.system("mkdir -p %s" % sample_set)
        os.chdir(sample_set)

        mutations = CollectionVCF(vcf_file=sample_dir + sample_set + suffix,
                                  from_file=True)
        print("Totaly %s mutations" % len(mutations))

        mutations.test_thresholds(extracting_method='distance', threshold=(50, 5000, 100),
                                  testing_dir="testing_threshold")