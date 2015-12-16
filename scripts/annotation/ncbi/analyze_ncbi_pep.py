#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse

from CustomCollections.GeneralCollections import SynDict
parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with ncbi peptides  ")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files. By default - species name, with spaces replaced by underline")
parser.add_argument("-s", "--species_name", action="store", dest="species_name", required=True,
                    help="Species name as listed in pep file. Example: 'Felis catus'")

args = parser.parse_args()

if not args.output_prefix:
    args.output_prefix = args.species_name.replace(" ", "_")

pep_description_file = "%s.pep.description" % args.output_prefix
pep_uniq_description_file = "%s.pep.description.uniq" % args.output_prefix
pep_uniq_description_no_isoform_versions = "%s.pep.sorted.description.no_isoform_versions" % args.output_prefix
pep_description_collapsed_isoforms = "%s.pep.collapsed_isoforms.description" % args.output_prefix

get_pep_decription_str = "grep -P '^>' %s | sed 's/^>//;s/\[%s\]//;s/ /\\t/' | sort -t $'\\t' -k 2 -k 1 > %s" % (args.input,
                                                                                                                  args.species_name,
                                                                                                                  pep_description_file)
get_uniq_descriptions_str = "uniq -f 1 %s > %s" % (pep_description_file, pep_uniq_description_file)
remove_isoform_versions_str = "sed s/isoform.*// %s > %s" % (pep_uniq_description_file,
                                                             pep_uniq_description_no_isoform_versions)

for exe_string in get_pep_decription_str, get_uniq_descriptions_str, remove_isoform_versions_str:
    print(exe_string)
    os.system(exe_string)


syn_dict = SynDict()
syn_dict.read(pep_uniq_description_no_isoform_versions, header=False, separator="\t", allow_repeats_of_key=True,
              split_values=True, values_separator=",", key_index=0, value_index=1, comments_prefix="#")
syn_dict.write(pep_description_collapsed_isoforms, splited_values=True, values_separator=",")