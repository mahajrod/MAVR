#!/usr/bin/env python
import os
import re
import time
import xmltodict

from collections import Iterable

import numpy as np

from collections import OrderedDict
from Bio import SeqIO, Entrez
from Bio.SeqRecord import SeqRecord
from Routines import FileRoutines
from CustomCollections.GeneralCollections import IdList, SynDict, TwoLvlDict, IdSet

from urllib2 import URLError


class EnsemblRoutines(FileRoutines):
    def __init__(self):
        FileRoutines.__init__(self)

    @staticmethod
    def convert_biomart_protein_annotation_to_gff(biomart_annotations, output_gff, separator="\t",
                                                  extraction_mode="pfam"):

        in_fd = open(biomart_annotations, "r")
        protein_id_header = "Protein ID"

        if extraction_mode == "pfam":
            domain_id_header = "Pfam ID"
            domain_start_header = "Pfam start"
            domain_end_header = "Pfam end"

        header_list = in_fd.readline().strip().split(separator)

        for header_element in  protein_id_header, domain_id_header, domain_start_header, domain_end_header:
            if header_element not in header_list:
                raise ValueError("'%s' is absent in header of the file" % header_element)

        protein_id_index = header_list.index(protein_id_header)
        domain_id_index = header_list.index(domain_id_header)
        domain_start_index = header_list.index(domain_start_header)
        domain_end_index = header_list.index(domain_end_header)

        with open(output_gff, "w") as out_fd:
            for line in in_fd:
                tmp_list = line.strip().split(separator)
                protein_id = tmp_list[protein_id_index]
                domain_id = tmp_list[domain_id_index]
                domain_start = tmp_list[domain_start_index]
                domain_end = tmp_list[domain_end_index]

                out_fd.write("%s\tensembl\tdomain\t%s\t%s\t.\t+\t.\tID=%s; Description=%s\t" % (protein_id,
                                                                                                domain_start,
                                                                                                domain_end,
                                                                                                domain_id,
                                                                                                domain_id))

        in_fd.close()