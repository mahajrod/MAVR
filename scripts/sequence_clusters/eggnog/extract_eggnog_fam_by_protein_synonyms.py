#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import EggNOGRoutines
from RouToolPa.Collections.General import SynDict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input file with EggNOG groups")
parser.add_argument("-y", "--protein_common_names_synonyms", action="store", dest="protein_common_names_synonyms",
                    required=True,
                    help="File with common names and associated ensembl ids. "
                         "Format of lines: common name\t<list of ensembl ids>")
parser.add_argument("-s", "--species_id", action="store", dest="species_id",
                    help="NCBI taxonomy id of the species. "
                         "If not set ensembl ids will be treated as already labeled by species id")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Output_prefix")
parser.add_argument("-r", "--allow_repeating_of_keys", action="store_true", dest="allow_repeating_of_keys",
                    help="Allow repeating of keys in synonym file")

args = parser.parse_args()

EggNOGRoutines.extract_eggnog_fam_by_protein_syn_dict(SynDict(filename=args.input, split_values=True),
                                                      SynDict(filename=args.protein_common_names_synonyms,
                                                              split_values=True,
                                                              allow_repeats_of_key=args.allow_repeating_of_keys),
                                                      output_prefix=args.output_prefix,
                                                      species_id=args.species_id)
