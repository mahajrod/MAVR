#!/usr/bin/env python2

__author__ = 'mahajrod'

import os
from Parsers.Cufflinks import CollectionFPKMTracking
from Parsers.GFF import CollectionGFF

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/combined_vcf/"

    #bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all.gff"
    bad_regions_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/masked_regions/LAN210_v0.10m_masked_all_not_in_good_genes.gff"
    bad_regions = CollectionGFF(input_file=bad_regions_file,
                                from_file=True)

    """
    expression_data_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/expression/Yeast_RNA_Seq/"
    expression_samples_dir_list = ["S1_Nagal_clout/",
                                   "S1_Yassour_clout/",
                                   "S2_Nagal_clout/",
                                   "S2_Yassour_clout/"]
    """
    expression_data_dir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/expression/stranded_expression/"
    expression_samples_dir_list = ["S1_Nagal_clout/",
                                   "S1_Yassour_clout/",
                                   "S2_Nagal_clout/",
                                   "S2_Yassour_clout/",
                                   "S3_Nagal_clout/",
                                   "S3_Yassour_clout/",
                                   "S4_Nagal_clout/",
                                   "S4_Yassour_clout/",
                                   "S5_Nagal_clout/",
                                   "S5_Yassour_clout/",
                                   "S6_Nagal_clout/",
                                   "S6_Yassour_clout/"]



    for expression_sample_dir in expression_samples_dir_list:
        os.chdir(expression_data_dir)
        os.chdir(expression_sample_dir)
        gene_expresion_data = CollectionFPKMTracking(from_file=True, input_file="genes.fpkm_tracking")

        expression = "((bad_region.start <= self.pos <= bad_region.end) \
                           or (bad_region.start <= self.pos + self.length - 1 <= bad_region.end) \
                           or (self.pos <= bad_region.start <= self.pos + self.length - 1))"

        print("Handling %s" % expression_sample_dir[:-1])
        gene_expresion_data.check_location(bad_regions, expression=expression)
        filtered_records, filtered_out_records = gene_expresion_data.filter_by_expression("\'BR\' not in record.flags")
        filtered_records.write("filtered_genes.fpkm_tracking")
        filtered_out_records.write("filtered_out_genes.fpkm_tracking")