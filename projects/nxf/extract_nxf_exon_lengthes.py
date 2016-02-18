#!/usr/bin/env python
import os
from Bio import SeqIO

from Routines import FileRoutines

workdir = "/home/mahajrod/Genetics/Projects/nxf/nxf_arthropoda/"
data_dir = "/home/mahajrod/Genetics/Projects/nxf/nxf_arthropoda/data/"

os.chdir(workdir)

data_files = FileRoutines.make_list_of_path_to_files([data_dir])

record_dict = SeqIO.index_db("tmp.idx", data_files, format="genbank")

for record_id in record_dict:
    for feature in record_dict[record_id].features:
        if feature.type == "mRNA":
            print("AAAAAAAAAAA")
            print(",".join(record_dict[record_id].annotations["taxonomy"]))
            print(record_dict[record_id].annotations["organism"])
            print(feature.qualifiers["product"])
            print feature.sub_features
            print feature.location
            print type(feature.location)
            #print feature.location.start
            for location in feature.location.parts:
                print location

