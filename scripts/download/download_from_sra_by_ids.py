#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from Bio import SeqIO
from Tools.Annotation import SNPeff
from Routines.File import make_list_of_path_to_files


def path_from_id(entry_id):
    ncbi_ftp = "ftp://ftp-trace.ncbi.nlm.nih.gov/"
    sra_reads_dir = "sra/sra-instant/reads/"
    id_code_dict = {
                "DRR": "ByRun",
                "ERR": "ByRun",
                "SRR": "ByRun",

                "DRX": "ByExp",
                "ERX": "ByExp",
                "SRX": "ByExp",

                "DRS": "BySample",
                "ERS": "BySample",
                "SRS": "BySample",

                "DRP": "ByStudy",
                "ERP": "ByStudy",
                "SRP": "ByStudy"
                }
    id_group = entry_id[:3]
    id_subgroup = entry_id[:6]

    id_type = id_code_dict[id_group]

    return "%s%s%s/sra/%s/%s/%s/" % (ncbi_ftp, sra_reads_dir, id_type, id_group, id_subgroup, entry_id)

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--ids", action="store", dest="ids", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of SRA ids to download")

args = parser.parse_args()

for entry_id in args.ids:
    ftp_path = path_from_id(entry_id)
    os.system("wget --no-host-directories -rc -t 500 %s" % ftp_path)
