#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from RouToolPa.Collections.General import TwoLvlDict, IdList, IdSet, SynDict
from RouToolPa.Routines.File import check_path
from RouToolPa.Routines.File import read_synonyms_dict




def filter_nonassembled(families):
    for entry in set(families):
        if entry[0] == "C" or entry[0] == "I" or entry[0] == "M":
            return False
    return True


def filter_splited_to_several_fam(families):
    for entry in set(families):
        if ";" in entry:
            return False
        if entry[0] != "I" and "_" in entry:
            return False
    return True

def filter_different_assembly(families):
    return False if len(set(families)) > 1 else True


parser = argparse.ArgumentParser()

parser.add_argument("-s", "--species_list", action="store", dest="species_list", type=lambda s: s.split(","),
                    required=True,
                    help="Comma-separated list of species")
parser.add_argument("-d", "--species_dir", action="store", dest="species_dir", default="./", type=check_path,
                    help="Comma-separated list of species")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file. Default: stdout")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

species_syn_dict = TwoLvlDict()

for species in args.species_list:
    species_syn_dict[species] = read_synonyms_dict("%s%s/all.t" % (args.species_dir, species))

species_syn_dict.write("families_all_species.t", absent_symbol=".")

nonassembled = species_syn_dict.filter_by_line(filter_nonassembled)
species_syn_dict.write("correctly_assembled_families_species.t", absent_symbol=".")


nonassembled.write("not_assembled_families_in_all_species.t", absent_symbol=".")
complicated_families_dict = nonassembled.filter_by_line(filter_splited_to_several_fam)
complicated_families_dict.write("complicated_families.t", absent_symbol=".")

complicated_families_syn_dict = SynDict()
complicated_families_syn_ids = IdSet()
sl_keys = list(complicated_families_dict.sl_keys())
for sl_key in sl_keys:
    sp_set = set()
    for species in complicated_families_dict:
        if sl_key not in complicated_families_dict[species]:
            continue
        tmp = complicated_families_dict[species][sl_key].split(";")
        for i in range(0, len(tmp)):
            if "_" in tmp[i]:
                tmp[i] = tmp[i][2:]
            tmp[i] = tmp[i].split(",")
            for syn_id in tmp[i]:
                complicated_families_syn_ids.add(syn_id)
                sp_set.add(syn_id)
    complicated_families_syn_dict[sl_key] = sp_set
complicated_families_syn_dict.write("complicated_families_connections.t", splited_values=True)

for entry in complicated_families_dict.all_values():
    tmp = entry.split(";")
    for i in range(0, len(tmp)):
        if "_" in tmp[i]:
            tmp[i] = tmp[i][2]
        tmp[i] = tmp[i].split(",")
        for syn_id in tmp[i]:
            complicated_families_syn_ids.add(syn_id)
complicated_families_syn_ids.write("complicated_families_check.ids")

nonassembled.write("splited_to_several_families.t", absent_symbol=".")

assemled_to_different_families = species_syn_dict.filter_by_line(filter_different_assembly)
species_syn_dict.write("correctly_assembled_families_in_all_species.t", absent_symbol=".")
assemled_to_different_families.write("assembled_to_different_families_in_all_species.t", absent_symbol=".")

correctly_assembled_families_synonym = IdList(set(species_syn_dict.all_values()))
assemled_to_different_families_synonym = IdList(set(assemled_to_different_families.all_values()))

correctly_assembled_families_synonym.write("correctly_assembled_families_syn_in_all_species.ids")
assemled_to_different_families_synonym.write("assembled_to_different_families_syn_in_all_species.ids")
if args.output != "output":
    out_fd.close()
