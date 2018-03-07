#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from collections import OrderedDict

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input vcf file")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

vep_header = "#Chrom\tpos\tid\tref\talt\tqual\tfilter\tEAS_AF\tEUR_AF\tAFR_AF\tAMR_AF\tSAS_AF\t" \
         "Allele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\t" \
         "CDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\t" \
         "DISTANCE\tSTRAND\tSIFT\tPolyPhen\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\n"

erb_header = "#Chrom\tpos\tid\tref\talt\tqual\tfilter\tEAS_AF\tEUR_AF\tAFR_AF\tAMR_AF\tSAS_AF\t" \
            "Allele\t_Gene\t_Feature\t_Feature_type\t_Consequence\n"

with open(args.input, "r") as in_fd:
    with open("%s.vep.tab" % args.output_prefix, "w") as out_fd:
        with open("%s.erb.tab" % args.output_prefix, "w") as erb_fd:
            out_fd.write(vep_header)
            erb_fd.write(erb_header)
            for line in in_fd:
                if line[0] == "#":
                    continue

                line_list = line.strip().split("\t")
                info_dict = OrderedDict()
                for entry in map(lambda e: e.split("="), line_list[7].split(";")):
                    #print entry
                    if len(entry) == 2:
                        info_dict[entry[0]] = entry[1]
                #    info_dict = OrderedDict(map(lambda e: e.split("="), line_list[7].split(";")))
                for key in "CSQ", "ERB":
                    if key in info_dict:
                        info_dict[key] = map(lambda e: e.split("|"), info_dict[key].split(","))
                if "CSQ" not in info_dict:
                    continue

                out_string = "%s\t%s\t%s\t%s\t%s\t%s\t" % ("\t".join(line_list[0:7],),
                                                           info_dict["EAS_AF"],
                                                           info_dict["EUR_AF"],
                                                           info_dict["AFR_AF"],
                                                           info_dict["AMR_AF"],
                                                           info_dict["SAS_AF"])
                for entry in info_dict["CSQ"]:
                    for element in entry:
                        if element != "":
                            break
                    else:
                        continue

                    full_string = out_string + "\t".join(entry) + "\n"
                    out_fd.write(full_string)
                if "ERB" in info_dict:
                    for entry in info_dict["ERB"]:
                        for element in entry:
                            if element != "":
                                break
                        else:
                            continue

                        full_string = out_string + "\t".join(entry) + "\n"
                        erb_fd.write(full_string)

