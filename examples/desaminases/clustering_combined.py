#!/usr/bin/env python2
import os

from BCBio import GFF

from Parsers.VCF import ReferenceGenome, CollectionVCF, ref_alt_variants
from Parsers.CCF import CollectionCCF
from Parsers.GFF import CollectionGFF
if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.idx")

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
    """
    sample_set_names_list = [#"HAP",
                             "HAP_sub1"]
    """


    clustering_dir = "clustering"
    rainfall_dir = "rainfall"
    distance_threshold = 1000
    reference.find_gaps()
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))
    min_cluster_size = 3
    #bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all.gff"

    bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all_not_in_good_genes.gff"
    bad_regions = CollectionGFF(input_file=bad_regions_file,
                                from_file=True)
    gff_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    annotations_dict = {}
    annotation_synonym_dict = {"three_prime_UTR": "3'_UTR",
                               "five_prime_UTR": "5'_UTR",
                               "snoRNA": "ncRNA",
                               "snRNA": "ncRNA"
                               }
    with open(gff_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            annotations_dict[record.id] = record

    bad_region_dict = {}
    with open(bad_regions_file) as gff_fd:
        for record in GFF.parse(gff_fd):
            bad_region_dict[record.id] = record

    for sample_set_name in sample_set_names_list:
        print("Handling %s" % sample_set_name)

        os.chdir(workdir)
        os.system("mkdir -p %s" % sample_set_name)
        os.chdir(sample_set_name)
        os.system("mkdir -p %s %s" % (clustering_dir, rainfall_dir))
        #os.system("pwd")
        mutations = CollectionVCF(vcf_file="../SNP_annotated_raw_vcf/%s_SNP.vcf" % sample_set_name,
                                  from_file=True)

        """
        mutations.rainfall_plot("%s_mutations" % (sample_set_name), ref_genome=reference, draw_gaps=True,
                                masked_regions=bad_region_dict)
        """

        mutations.get_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict)
        mutations.check_location(bad_regions)
        mutations.check_by_ref_and_alt(ref_alt_variants["desaminases"], "DA")



        annotation_black_list = ["gene", "region", "ARS", "long_terminal_repeat",
                                 "noncoding_exon", "intron", "repeat_region", "telomere", "gene_cassette",
                                 "five_prime_UTR_intron"]
        """
        mutations.location_pie(annotation_black_list=annotation_black_list,
                               figsize=(30, 30),
                               pie_filename="%s_SNP_raw_variant_location_pie.svg" % sample_set_name,
                               full_genome_pie_filename="%s_SNP_raw_variant_location_full_genome_pie.svg" % sample_set_name,
                               counts_filename="%s_SNP_raw_variant_location_counts.t" % sample_set_name)
        """
        print("Totaly %s mutations" % len(mutations))

        #mutations.hierarchical_clustering(sample_name=sample, save=True, clustering_dir=clustering_dir)
        sample_set_name_adjusted = sample_set_name + "_adjusted"
        clusters = mutations.get_clusters(sample_name=sample_set_name, save_clustering=True,
                                          extracting_method="distance",
                                          threshold=distance_threshold, cluster_distance='average',
                                          dendrogramm_max_y=2000, dendrogramm_color_threshold=1000,
                                          clustering_dir=clustering_dir, split_by_regions=False)
        clusters.subclustering()
        clusters.statistics(filename="%s/%s_cluster_size_distribution.svg" % (clustering_dir, sample_set_name))
        clusters.adjust(border_limit=None, min_size_to_adjust=1, remove_border_subclusters=True, remove_size_limit=1)
        clusters.statistics(filename="%s/%s_cluster_size_distribution.svg" % (clustering_dir, sample_set_name_adjusted))
        clusters.check_location(bad_regions)
        clusters.write("%s/%s_raw.ccf" % (clustering_dir, sample_set_name))

        clusters.check_flags(["DA"], mismatch_list=[1], expression_list=["record.count_samples() <= 1"],
                             remove_mismatch_list=[True])
        clusters.get_location(annotations_dict, use_synonym=True, synonym_dict=annotation_synonym_dict)
        clusters.write("%s/%s.ccf" % (clustering_dir, sample_set_name_adjusted))

        filtered_clusters, filtered_out_clusters = clusters.filter_by_flags(black_flag_list=["IP", "BR"])
        filtered_clusters.write("%s/%s_not_in_br_no_id.ccf" % (clustering_dir, sample_set_name_adjusted))
        """
        filtered_clusters.location_pie(annotation_black_list=annotation_black_list,
                                       pie_filename="%s_good_cluster_location_pie.svg" % sample_set_name_adjusted,
                                       full_genome_pie_filename="%s_good_cluster_full_genome_pie.svg" % sample_set_name_adjusted,
                                       counts_filename="%s_good_cluster_location_counts.t" % sample_set_name_adjusted)
        """
        filtered_out_clusters.write("%s/%s_in_br_id.ccf" % (clustering_dir, sample_set_name_adjusted))
        filtered_clusters.statistics(filename="%s/%s_not_in_br_no_id_cluster_size_distribution.svg" % (clustering_dir, sample_set_name_adjusted))
        filtered_out_clusters.statistics(filename="%s/%s_in_br_id_cluster_size_distribution.svg" % (clustering_dir, sample_set_name_adjusted))

        if "HAP" not in sample_set_name:
            filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_flags(white_flag_list=["DA"])
            filtered_clusters.write("%s/%s_not_in_br_no_id_da.ccf" % (clustering_dir, sample_set_name_adjusted))
            """
            filtered_clusters.location_pie(annotation_black_list=annotation_black_list,
                                           pie_filename="%s_good_cluster_desaminase_location_pie.svg" % sample_set_name_adjusted,
                                           full_genome_pie_filename="%s_good_cluster_desaminase_full_genome_pie.svg" % sample_set_name_adjusted,
                                           counts_filename="%s_good_cluster_desaminase_location_counts.t" % sample_set_name_adjusted)
            """
            filtered_out_clusters.write("%s/%s_not_in_br_no_id_non_da.ccf" % (clustering_dir, sample_set_name_adjusted))
            filtered_clusters.check_strandness()

            filtered_clusters.statistics(filename="%s/%s_not_in_br_no_id_da_cluster_size_distribution.svg" % (clustering_dir,sample_set_name_adjusted))
            filtered_out_clusters.statistics(filename="%s/%s_not_in_br_no_id_non_da_distribution.svg" % (clustering_dir, sample_set_name_adjusted))

            filtered_clusters.heatmap_statistics(filename="%s/%s_not_in_br_no_id_da_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power", "Homogeneity"))
            #filtered_clusters.strandness_statistics(filename="%s/%s_not_in_br_no_id_da_strandness_statistics.svg" % (clustering_dir,sample_set_name_adjusted))
        else:
            filtered_clusters.heatmap_statistics(filename="%s/%s_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power"))


        cluster_mutations = filtered_clusters.extract_vcf()
        cluster_mutations.write("%s/%s_cluster_mutations.vcf" % (clustering_dir, sample_set_name_adjusted))
        """
        cluster_mutations.rainfall_plot("%s_cluster_mutations" % (sample_set_name_adjusted), ref_genome=reference, draw_gaps=True,
                                        masked_regions=bad_region_dict)

        cluster_mutations.location_pie(annotation_black_list=annotation_black_list,
                                       figsize=(30, 30),
                                       pie_filename="%s_good_cluster_mutations_location_pie.svg" % sample_set_name,
                                        full_genome_pie_filename="%s_good_cluster_mutations_location_full_genome_pie.svg" % sample_set_name,
                                       counts_filename="%s_good_cluster_mutations_location_counts.t" % sample_set_name)
        """

        filtered_clusters, filtered_out_clusters = filtered_clusters.filter_by_size(min_size=min_cluster_size)
        filtered_clusters.write("%s/%s_size_3+.ccf" % (clustering_dir, sample_set_name_adjusted))
        """
        filtered_clusters.location_pie(annotation_black_list=annotation_black_list,
                                       pie_filename="%s_3+_cluster_location_pie.svg" % sample_set_name_adjusted,
                                       counts_filename="%s_3+_cluster_location_counts.t" % sample_set_name_adjusted,
                                       full_genome_pie_filename="%s_3+_cluster_location_full_genome_pie.svg" % sample_set_name_adjusted,)
        """
        filtered_out_clusters.write("%s/%s_size_less_3.ccf" % (clustering_dir, sample_set_name_adjusted))

        if "HAP" not in sample_set_name:
            #filtered_clusters.strandness_statistics(filename="%s/%s_3+_not_in_br_no_id_da_strandness_statistics.svg" % (clustering_dir,sample_set_name_adjusted))
            #filtered_out_clusters.strandness_statistics(filename="%s/%s_less_3_not_in_br_no_id_non_da_strandness_statistics.svg" % (clustering_dir, sample_set_name_adjusted))
            filtered_clusters.heatmap_statistics(filename="%s/%s_3+_not_in_br_no_id_da_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power", "Homogeneity"))
            filtered_out_clusters.heatmap_statistics(filename="%s/%s_less_3_not_in_br_no_id_non_da_heatmap_statistics.svg" % (clustering_dir, sample_set_name_adjusted),
                                                     additional_data=("Median", "Mean", "Power", "Homogeneity"))
        else:
            filtered_clusters.heatmap_statistics(filename="%s/%s_3+_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir,sample_set_name_adjusted),
                                                 additional_data=("Median", "Mean", "Power"))
            filtered_out_clusters.heatmap_statistics(filename="%s/%s_less_3_not_in_br_no_id_heatmap_statistics.svg" % (clustering_dir, sample_set_name_adjusted),
                                                     additional_data=("Median", "Mean", "Power"))

        cluster_mutations = filtered_clusters.extract_vcf()
        cluster_mutations.write("%s/%s_3+_cluster_mutations.vcf" % (clustering_dir, sample_set_name_adjusted))
        """
        cluster_mutations.rainfall_plot("%s_3+_cluster_mutations" % (sample_set_name_adjusted), ref_genome=reference, draw_gaps=True,
                                        masked_regions=bad_region_dict)

        cluster_mutations.location_pie(annotation_black_list=annotation_black_list,
                                       figsize=(30, 30),
                                       pie_filename="%s_3+_good_cluster_mutations_location_pie.svg" % sample_set_name,
                                       full_genome_pie_filename="%s_3+_good_cluster_mutations_location_full_genome_pie.svg" % sample_set_name,
                                       counts_filename="%s_3+good_cluster_mutations_location_counts.t" % sample_set_name)
        """