import os

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/clusters/all/all/"

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
    annotations = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/annotations/merged_annotations_Nagalakshmi_tranf_to_LAN210_v0.10m.gff3"
    os.chdir(workdir)
    right = 50
    left = 200
    for sample_set in sample_set_names_list:
        print("Handling %s" % sample_set)
        vcf_file = "%s_good.vcf" % sample_set
        start_hist_prefix = "%s_start_hist_r_%i_l_%i" % (sample_set, right, left)
        end_hist_prefix = "%s_end_hist_r_%i_l_%i" % (sample_set, right, left)
        gene_variants = "%s_gene_variants_r_%i_l_%i.t" % (sample_set, right, left)
        os.system("CDS_region_variants.py -f %s -v %s -s %s -e %s -l %i -r %i -g %s" %
                  (annotations, vcf_file, start_hist_prefix, end_hist_prefix, left, right, gene_variants))