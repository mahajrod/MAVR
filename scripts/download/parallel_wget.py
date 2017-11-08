#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys
import argparse

from CustomCollections.GeneralCollections import IdList

from multiprocessing import Pool

#from Tools.Abstract import Tool
#from Routines import NCBIRoutines

from Tools.LinuxTools import Wget

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--link_file", action="store", dest="link_file",
                    help="File with links")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of simultaneous downloads")

args = parser.parse_args()

loader = IdList()
link_list = loader.read(args.link_file)

Wget.threads = args.threads
Wget.parallel_download(link_list)

"""
options_list = []
for entry_id in id_list:
    ftp_path = NCBIRoutines.get_sra_ftp_path_from_id(entry_id)
    options_list.append("--no-host-directories -rc -t 500 %s" % ftp_path)


tool = Tool(cmd="wget", max_threads=args.threads)

tool.parallel_execute(options_list)
"""