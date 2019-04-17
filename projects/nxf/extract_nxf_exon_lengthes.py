#!/usr/bin/env python
import os
from Bio import SeqIO
from RouToolPa.Routines import FileRoutines


workdir = "/home/mahajrod/Genetics/Projects/nxf/nxf_arthropoda/"
data_dir = "/home/mahajrod/Genetics/Projects/nxf/nxf_arthropoda/data/"

os.chdir(workdir)

data_files = FileRoutines.make_list_of_path_to_files([data_dir])

record_dict = SeqIO.index_db("tmp.idx", data_files, format="genbank")

print("#organism\ttaxonomy\tregion_id\ttranscript_id\tproduct\texon_len")
for record_id in record_dict:
    for feature in record_dict[record_id].features:
        if feature.type == "mRNA":
            mRNA_string = ""
            mRNA_string += "%s" % record_dict[record_id].annotations["organism"]
            mRNA_string += "\t%s" % (";".join(record_dict[record_id].annotations["taxonomy"]))
            mRNA_string += "\t%s" % record_id
            mRNA_string += "\t%s" % (feature.qualifiers["transcript_id"][0] if "transcript_id" in feature.qualifiers else ".")
            mRNA_string += "\t%s" % (feature.qualifiers["product"][0] if "product" in feature.qualifiers else ".")

            location_lenths = []

            for location in feature.location.parts:
                location_lenths.append(len(location))

            mRNA_string += "\t%s" % (",".join(map(str, location_lenths)))

            print(mRNA_string)

            #print(feature.qualifiers)
