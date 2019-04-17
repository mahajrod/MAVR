#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Expression import Subread

parser = argparse.ArgumentParser()

parser.add_argument("-a", "--alignment", action="store", dest="alignment", required=True,
                    help="Alignment file")
parser.add_argument("-e", "--gff_for_subread", action="store", required=True,
                    dest="gff_for_subread", help="Gff file with annotations for Subread")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    required=True, help="Prefix of output files")
parser.add_argument("-m", "--min_read_fraction_overlap", action="store", default=1.0, type=float,
                    dest="min_read_fraction_overlap",
                    help="Default read fraction overlap. Default - 1.0(100%%)")
parser.add_argument("--feature_type", action="store", default="exon",
                    dest="feature_type",
                    help="Feature type from annotation gff to be used for counting reads. Default - 'exon'")
parser.add_argument("--feature_id_attribute", action="store", default="gene_id",
                    dest="feature_id_attribute",
                    help="Feature id attribute from annotation gff to be considered for counting reads. "
                         "Default - 'gene_id'")
parser.add_argument("-s", "--sample", action="store", dest="sample",
                    help="Sample name. By default filename from subread file is used")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use")

args = parser.parse_args()

Subread.threads = args.threads

Subread.count_miRNA_reads(args.alignment, args.gff_for_subread, args.output_prefix, annotation_file_type="GTF",
                          min_read_fraction_overlap=args.min_read_fraction_overlap,
                          feature_type_to_use=args.feature_type,
                          attribute_type_to_use=args.feature_id_attribute,
                          sample_name=args.sample)
