#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from CustomCollections.GeneralCollections import IdList

from multiprocessing import Pool

from Tools.Abstract import Tool


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

parser.add_argument("-i", "--ids", action="store", dest="ids",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of SRA ids to download")
parser.add_argument("-f", "--id_file", action="store", dest="id_file",
                    help="File with SRA ids(one per line) to download")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of simultaneous downloads")

args = parser.parse_args()

if (not args.ids) and (not args.id_file):
    raise ValueError("Both ids and id file were not set")

loader = IdList()
id_list = loader.read(args.id_file) if args.id_file else args.ids

options_list = []
for entry_id in id_list:
    ftp_path = path_from_id(entry_id)
    options_list.append("--no-host-directories -rc -t 500 %s" % ftp_path)


tool = Tool(cmd="wget", max_threads=args.threads)

tool.parallel_execute(options_list)
