#!/usr/bin/env python
import os

from BCBio import GFF

from Parser.VCF import ReferenceGenome, CollectionVCF
from Parser.CCF import CollectionCCF
from Parser.GFF import CollectionGFF
if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all"

    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.idx")

    suffix = "_GATK_best_merged.vcf"
    clustering_dir = "clustering"
    alignment_dir = "alignment_LAN210_v0.10m"
    distance_threshold = 1000
    reference.find_gaps()
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))
    min_cluster_size = 3

    bad_regions = CollectionGFF(input_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/masked_regions/LAN210_v0.9m_masked_all.gff",
                                from_file=True)
    gff_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    annotations_dict = {}
    with open(gff_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            annotations_dict[record.id] = record

    for sample in samples_list:
        print("Handling %s" % sample)

        os.chdir(workdir)
        os.chdir(sample)
        if alignment_dir not in os.listdir("."):
            continue
        os.chdir(alignment_dir)
        os.system("mkdir -p %s" % clustering_dir)
        mutations = CollectionVCF(vcf_file=sample + suffix,
                                  from_file=True)
        mutations.get_location(annotations_dict)
        #for record in mutations:
        #    print(record.description)
        mutations.location_pie(annotation_black_list=["gene"], figsize=(40, 40),
                               pie_filename="variant_location_pie.svg",
                               counts_filename="variant_location_counts.t")
        print("Totaly %s mutations" % len(mutations))

        #mutations.hierarchical_clustering(sample_name=sample, save=True, clustering_dir=clustering_dir)

        clusters = mutations.get_clusters(sample_name=sample, save_clustering=True, extracting_method="distance",
                                          threshold=distance_threshold, cluster_distance='average',
                                          dendrogramm_max_y=2000, dendrogramm_color_threshold=1000,
                                          clustering_dir=clustering_dir, split_by_regions=False)
        clusters.check_record_location(bad_regions)
        clusters.get_location(annotations_dict)
        clusters.write("%s/%s%s.ccf" % (clustering_dir, sample, suffix))

        filtered_clusters, filtered_out_clusters = clusters.filter_by_size(min_size=min_cluster_size)
        filtered_clusters.write("%s/%s%s_size_3+.ccf" % (clustering_dir, sample, suffix))
        filtered_clusters.location_pie(annotation_black_list=["gene"],
                                       pie_filename="3+_cluster_location_pie.svg",
                                       counts_filename="3+_cluster_location_counts.t")
        filtered_out_clusters.write("%s/%s%s_size_less_3.ccf" % (clustering_dir, sample, suffix))

        filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_flags(black_flag_list=["IP", "BR"])
        filtered_clusters.write("%s/%s%s_size_3+_not_in_br_no_id.ccf" % (clustering_dir, sample, suffix))
        filtered_clusters.location_pie(annotation_black_list=["gene"],
                                       pie_filename="3+_good_cluster_location_pie.svg",
                                       counts_filename="3+_good__cluster_location_counts.t")
        filtered_out_clusters.write("%s/%s%s_size_3+_in_br_id.ccf" % (clustering_dir, sample, suffix))
        #mutations.test_thresholds(save_clustering=True, testing_dir="testing_theshold_inconsistent")
        #mutations.test_thresholds(extracting_method='distance', threshold=(50, 5000, 100),
        #                          testing_dir="testing_threshold")
