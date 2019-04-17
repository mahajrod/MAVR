#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Collections.General import SynDict, IdSet


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input", required=True,
                    help="Input file with pep families")
parser.add_argument("-o", "--output_file", action="store", dest="output", required=True,
                    help="Output cds fam file")
#parser.add_argument("-d", "--output_dir", action="store", dest="output_dir", required=True,
#                    help="Output directory with cds for")
#parser.add_argument("-c", "--cds_dir", action="store", dest="cds_dir", required=True,
#                    help="Directory with cds files. Filenames must look like <species>.cds")
parser.add_argument("-a", "--accordance_dir", action="store", dest="accordance_dir", required=True,
                    help="Directory with cds-to-pep id accordance files. Filenames must look like <species>.accordance")
parser.add_argument("-f", "--families_with_errors", action="store", dest="fam_error", default="error.fam.ids",
                    help="File to write ids of families with errors")
parser.add_argument("-s", "--species_set", action="store", dest="species_set",
                    help="Comma-separated list of species.")

parser.add_argument("-l", "--name_last", action="store_false", dest="name_first", default=True,
                    help="Position of name of species in gene_id. Default: name first")
parser.add_argument("-e", "--name_separator", action="store", dest="name_separator", default=".",
                    help="Separator between species name and gene name. Default: '.'")
args = parser.parse_args()


args.species_set = set(args.species_set.split(","))

pep_fam_dict = SynDict()
pep_fam_dict.read(args.input, split_values=True)

cds_fam_dict = SynDict()

cds_dict = {}
accordance_dict = {}

for species in args.species_set:
    #cds_file = "%s/%s.cds" % (args.cds_dir, species)
    #cds_dict[species] = SeqIO.index(cds_file, format="fasta")

    accordance_file = "%s/%s.accordance" % (args.accordance_dir, species)
    accordance_dict[species] = SynDict()
    accordance_dict[species].read(accordance_file, key_index=1, value_index=0)


if args.name_first:
    def split_name(pep_name):
        gene_list = pep_name.split(args.name_separator)
        return gene_list[0], args.name_separator.join(gene_list[1:])
else:
    def split_name(pep_name):
        gene_list = pep_name.split(args.name_separator)
        return gene_list[-1], args.name_separator.join(gene_list[:-1])

families_with_errors = IdSet()
for family in pep_fam_dict:
    cds_fam_dict[family] = []
    for pep in pep_fam_dict[family]:
        species, pep_name = split_name(pep)
        if pep_name in accordance_dict[species]:
            cds_name = "%s%s%s" % (species, args.name_separator, accordance_dict[species][pep_name]) if args.name_first else \
                "%s%s%s" % (accordance_dict[species][pep_name], args.name_separator, species)
            cds_fam_dict[family].append(cds_name)
        else:
            print("%s %s %s doesn't have associated cds in accordance file" % (family, species, pep_name))
            families_with_errors.add(family)

for family in families_with_errors:
    cds_fam_dict.pop(family, None)

families_with_errors.write(args.fam_error)
cds_fam_dict.write(args.output, splited_values=True)
