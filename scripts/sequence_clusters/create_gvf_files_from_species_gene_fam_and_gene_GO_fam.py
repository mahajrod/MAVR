#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import SequenceClusterRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--species-gene_fam_file", action="store", dest="species_gene_fam_file", required=True,
                    help="Input species-gene fam file")
parser.add_argument("-g", "--gene-GO_fam_file", action="store", dest="gene_GO_fam_file", required=True,
                    help="Input gene-GO fam file")
#parser.add_argument("-c", "--cluster_column_index", action="store", dest="cluster_column_index", type=int, default=0,
#                    help="Index of cluster column in synonym file. Default: 0")
#parser.add_argument("-v", "--element_column_index", action="store", dest="element_column_index", type=int, default=1,
#                    help="Index of element column in synonym file.Default: 1")
#parser.add_argument("-e", "--separator", action="store", dest="column_separator", default='\t',
#                    help="Column separator in synonym file. Default: \\t")

parser.add_argument("-o", "--output_directory", action="store", dest="output_dir", required=True,
                    help="Output directory")

args = parser.parse_args()

SequenceClusterRoutines.create_gvf_files_from_species_gene_fam_and_gene_GO_fam(args.species_gene_fam_file,
                                                                               args.gene_GO_fam_file,
                                                                               args.output_dir)
