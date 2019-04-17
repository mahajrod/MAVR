#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Routines import SequenceRoutines
from RouToolPa.Collections.General import SynDict



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with ncbi peptides  ")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files. By default - species name, with spaces replaced by underline")
parser.add_argument("-s", "--species_name", action="store", dest="species_name", required=True,
                    help="Species name as listed in pep file. Example: 'Felis catus'")
parser.add_argument("-r", "--remove_predicted", action="store_true", dest="remove_predicted",
                    help="Remove predicted proteins")

args = parser.parse_args()

if not args.output_prefix:
    args.output_prefix = args.species_name.replace(" ", "_")

len_file = "%s.len" % args.input
pep_description_file = "%s.pep.description" % args.output_prefix
pep_uniq_description_file = "%s.pep.description.uniq" % args.output_prefix
pep_uniq_ids = "%s.pep.description.uniq.ids" % args.output_prefix
pep_uniq_description_no_isoform_versions = "%s.pep.sorted.description.no_isoform_versions" % args.output_prefix
pep_description_collapsed_isoforms = "%s.pep.collapsed_isoforms.description" % args.output_prefix
pep_description_collapsed_isoforms_with_len = "%s.pep.collapsed_isoforms.with_len.description" % args.output_prefix
pep_description_longest_isoform = "%s.pep.longest_isoform" % args.output_prefix
pep_description_longest_isoform_ids = "%s.pep.longest_isoform.ids" % args.output_prefix
pep_description_longest_isoform_pep = "%s.pep.longest_isoform.pep" % args.output_prefix
awk_extract_ids_string = "awk -F'\t' '{print $1}' %s > %s"

if args.remove_predicted:
    get_pep_decription_str = "grep -P '^>' %s | grep -v 'PREDICTED' | sed 's/^>//;s/\[%s\]//;s/ /\\t/' | sort -t '\t' -k 2 -k 1 > %s" % (args.input,
                                                                                                                                         args.species_name,
                                                                                                                                         pep_description_file)
else:
    get_pep_decription_str = "grep -P '^>' %s | sed 's/^>//;s/\[%s\]//;s/ /\\t/' | sort -t '\t' -k 2 -k 1 > %s" % (args.input,
                                                                                                                   args.species_name,
                                                                                                                   pep_description_file)
get_uniq_descriptions_str = "uniq -f 1 %s > %s" % (pep_description_file, pep_uniq_description_file)
remove_isoform_versions_str = "sed s/isoform.*// %s > %s" % (pep_uniq_description_file,
                                                             pep_uniq_description_no_isoform_versions)


for exe_string in get_pep_decription_str, get_uniq_descriptions_str, remove_isoform_versions_str:
    print(exe_string)
    os.system(exe_string)

os.system(awk_extract_ids_string % (pep_uniq_description_file, pep_uniq_ids))

syn_dict = SynDict()
syn_dict.read(pep_uniq_description_no_isoform_versions, header=False, separator="\t", allow_repeats_of_key=True,
              split_values=True, values_separator=",", key_index=1, value_index=0, comments_prefix="#")
syn_dict.write(pep_description_collapsed_isoforms, splited_values=True, values_separator=",")

length_dict = SequenceRoutines.get_lengths_from_seq_file(args.input, format="fasta", out_file=len_file)

descr_with_len_fd = open(pep_description_collapsed_isoforms_with_len, "w")
descr_longest_isoform_fd = open(pep_description_longest_isoform, "w")
descr_longest_isoform_ids_fd = open(pep_description_longest_isoform_ids, "w")

for gene in syn_dict:
    len_list = []
    longest_isoform = None
    max_len = 0
    for isoform_id in syn_dict[gene]:
        length = length_dict[isoform_id]
        len_list.append(length)
        if length > max_len:
            max_len = length
            longest_isoform = isoform_id

    descr_with_len_fd.write("%s\t%s\t%s\n" % (gene, ",".join(syn_dict[gene]), ",".join(map(str, len_list))))
    descr_longest_isoform_fd.write("%s\t%s\t%i\n" % (gene, longest_isoform, max_len))
    descr_longest_isoform_ids_fd.write(longest_isoform)
    descr_longest_isoform_ids_fd.write("\n")

for file_descriptor in descr_with_len_fd, descr_longest_isoform_fd, descr_longest_isoform_ids_fd:
    file_descriptor.close()

SequenceRoutines.extract_sequence_by_ids(args.input, pep_description_longest_isoform_ids,
                                         pep_description_longest_isoform_pep, format="fasta", verbose=True)
