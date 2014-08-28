#!/usr/bin/env python2
import os

from BCBio import GFF

from Parser.VCF import ReferenceGenome, CollectionVCF, ref_alt_variants
from Parser.CCF import CollectionCCF
from Parser.GFF import CollectionGFF
if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.idx")

    sample_set_names_list = ["PmCDA1_3d",
                         "PmCDA1_6d",
                         "A1_3d",
                         "A1_6d",
                         "A3G_3d",
                         "AID_3d",
                         "AID_6d",
                         "PmCDA1_sub1_3d",
                         "PmCDA1_sub1_6d"]

    clustering_dir = "clustering"
    distance_threshold = 1000
    reference.find_gaps()
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))
    min_cluster_size = 3

    bad_regions = CollectionGFF(input_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all.gff",
                                from_file=True)
    gff_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    annotations_dict = {}

    with open(gff_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            annotations_dict[record.id] = record

    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)

        os.chdir(workdir)
        os.system("mkdir -p %s" % sample_set_name)
        os.chdir(sample_set_name)
        os.system("mkdir -p %s" % clustering_dir)
        #os.system("pwd")
        mutations = CollectionVCF(vcf_file="../" + sample_set_name + ".vcf",
                                  from_file=True)
        mutations.get_location(annotations_dict)
        mutations.check_by_ref_and_alt(ref_alt_variants["desaminases"], "DA")
        #for record in mutations:
        #    print(record.description)
        annotation_black_list = ["gene", "region", "ARS", "long_terminal_repeat",
                                 "noncoding_exon", ]
        mutations.location_pie(annotation_black_list=annotation_black_list,
                               figsize=(30, 30),
                               pie_filename="variant_location_pie.svg",
                               full_genome_pie_filename="variant_location_full_genome_pie.svg",
                               counts_filename="variant_location_counts.t")
        print("Totaly %s mutations" % len(mutations))

        #mutations.hierarchical_clustering(sample_name=sample, save=True, clustering_dir=clustering_dir)

        clusters = mutations.get_clusters(sample_name=sample_set_name, save_clustering=True,
                                          extracting_method="distance",
                                          threshold=distance_threshold, cluster_distance='average',
                                          dendrogramm_max_y=2000, dendrogramm_color_threshold=1000,
                                          clustering_dir=clustering_dir, split_by_regions=False)
        clusters.check_record_location(bad_regions)
        clusters.get_location(annotations_dict)
        clusters.write("%s/%s.ccf" % (clustering_dir, sample_set_name))

        filtered_clusters, filtered_out_clusters = clusters.filter_by_size(min_size=min_cluster_size)
        filtered_clusters.write("%s/%s_size_3+.ccf" % (clustering_dir, sample_set_name))
        filtered_clusters.location_pie(annotation_black_list=annotation_black_list,
                                       pie_filename="3+_cluster_location_pie.svg",
                                       counts_filename="3+_cluster_location_counts.t",
                                       full_genome_pie_filename="3+_cluster_location_full_genome_pie.svg",)
        filtered_out_clusters.write("%s/%s_size_less_3.ccf" % (clustering_dir, sample_set_name))

        filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_flags(black_flag_list=["IP", "BR"])
        filtered_clusters.write("%s/%s_size_3+_not_in_br_no_id.ccf" % (clustering_dir, sample_set_name))
        filtered_clusters.location_pie(annotation_black_list=annotation_black_list,
                                       pie_filename="3+_good_cluster_location_pie.svg",
                                       full_genome_pie_filename="3+_good_cluster_full_genome_pie.svg",
                                       counts_filename="3+_good_cluster_location_counts.t")
        filtered_out_clusters.write("%s/%s_size_3+_in_br_id.ccf" % (clustering_dir, sample_set_name))

        filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_flags(white_flag_list=["DA"])
        filtered_clusters.write("%s/%s_size_3+_not_in_br_no_id_da.ccf" % (clustering_dir, sample_set_name))
        filtered_clusters.location_pie(annotation_black_list=annotation_black_list,
                                       pie_filename="3+_good_cluster_desaminase_location_pie.svg",
                                       full_genome_pie_filename="3+_good_cluster_desaminase_full_genome_pie.svg",
                                       counts_filename="3+_good_cluster_desaminase_location_counts.t")
        filtered_out_clusters.write("%s/%s_size_3+_nont_in_br_no_id_non_da.ccf" % (clustering_dir, sample_set_name))
