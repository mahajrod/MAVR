#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import sys
import argparse
from BCBio import GFF



parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_gff", action="store", dest="input_gff",
                    help="Gff file with annotations to extract")
parser.add_argument("-o", "--output_file", action="store", dest="output", default="stdout",
                    help="Output file with information about transcripts")
parser.add_argument("-l", "--longest_isoforms", action="store", dest="longest_isoforms",
                    help="File to write longest isoforms")

args = parser.parse_args()

out_fd = sys.stdout if args.output == "stdout" else open(args.output, "w")

for record in GFF.parse(open(args.input_gff)):
    for feature in record.features:
        #print feature
        if feature.type == "gene":
            transcript_id_list = []
            transcript_len_list = []
            CDS_len_list = []
            for subfeature in feature.sub_features:
                #print subfeature
                #print(subfeature.type)
                if subfeature.type == "mRNA" or subfeature.type == "transcript":
                    transcript_id_list.append(subfeature.id)
                    transcript_len_list.append(len(subfeature))
                    CDS_len = 0
                    for subsubfeature in subfeature.sub_features:
                        if subsubfeature.type == "CDS":
                            CDS_len += len(subsubfeature)
                    CDS_len_list.append(CDS_len)

            out_fd.write("%s\t%s\t%s\t%s\t%s\n" % (feature.id, ",".join(transcript_id_list),
                                                   ",".join(map(str, transcript_len_list)),
                                                   ",".join(map(str, CDS_len_list)),
                                                   ",".join(map(str, map(lambda x: x/3, CDS_len_list)))))

out_fd.close()



