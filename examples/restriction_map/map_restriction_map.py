#!/usr/bin/env python
import re, argparse, os
from collections import OrderedDict

from Bio import SeqIO


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_file", action="store", dest="input_file",
                    help="Input file with sequences")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix",
                    help="Prefix of output files")
parser.add_argument("-f", "--format", action="store", dest="format", default="fasta",
                    help="Format of input file. Default - fasta")

args = parser.parse_args()

sequence_dict = SeqIO.index_db("temp.idx", args.input_file, args.format)

unisexualis_map = "VCPXCSHSXPGCNHVPOPNDCONXXNBG"

restrictase_symbols_dict = {"B": "BamHI",
                            "C": "BclI",
                            "D": "SpeI",
                            "E": "EcoRI",
                            "G": "BglII",
                            "H": "HindIII",
                            "N": "NheI",
                            "O": "NcoI",
                            "P": "PvuII",
                            "S": "SacII",
                            "V": "EcoRV",
                            "X": "XbaI"
                            }

restrictase_names_dict = {"BamHI": "B",
                          "BclI": "C",
                           "SpeI": "D",
                           "EcoRI": "E",
                           "BglII": "G",
                           "HindIII": "H",
                           "NheI": "N",
                           "NcoI": "O",
                           "PvuII": "P",
                           "SacII": "S",
                           "EcoRV": "V",
                           "XbaI": "X"
                           }

restrictase_sites_dict = {
                    "BamHI": "GGATCC",
                    "BclI": "TGATCA",
                    "BglII": "AGATCT",
                    "EcoRI": "GAATTC",
                    "EcoRV": "GATATC",
                    "HindIII": "AAGCTT",
                    "NcoI": "CCATGG",
                    "NheI": "GCTAGC",
                    "PvuII": "CAGCTG",
                    "SacII": "CCGCGG",
                    "SpeI": "ACTAGT",
                    "XbaI": "TCTAGA"
                    }


def cart_sites(sequence, restrictase_sites_dict):
    sites_dictionary = {}
    for restrictase in restrictase_sites_dict:
        reg_exp = re.compile(restrictase_sites_dict[restrictase])
        mathes = reg_exp.finditer(sequence)
        for entry in mathes:
            if entry.start() in sites_dictionary:
                print("Warning: overlap of restriction sites")
            sites_dictionary[entry.start()] = restrictase_names_dict[restrictase]
    sorted_sites_dictionary = OrderedDict()
    sites_sequence = ""
    for position in sorted(sites_dictionary.keys()):
        sorted_sites_dictionary[position] = sites_dictionary[position]
        sites_sequence += sites_dictionary[position]
    return sorted_sites_dictionary, sites_sequence

sites_fd = open(args.output_prefix, "w")
sites_fd.write("#sequence\tsites\tsites_map\treverse_map\tsites_seq\n")
for record_id in sequence_dict:
    sites_dict, sites_string = cart_sites(str(sequence_dict[record_id].seq), restrictase_sites_dict)
    sites = "".join([str(position) + sites_dict[position] for position in sites_dict])
    sites_seq = ",".join([str(sequence_dict[record_id].seq[entry: entry+6]) for entry in sites_dict])
    sites_fd.write("%s\t%s\t%s\t%s\t%s\n" % (record_id, sites, sites_string if sites_string != "" else ".",
                                         sites_string[::-1] if sites_string != "" else ".", sites_seq))

sites_fd.close()
os.remove("temp.idx")