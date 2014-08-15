__author__ = 'mahajrod'

import os
from Parser.VCF import CollectionVCF, ReferenceGenome

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all"

    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.idx")
    #print(reference.reference_genome)
    reference.find_gaps()
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

        mutations.rainfall_plot(sample + suffix, ref_genome=reference, draw_gaps=True)