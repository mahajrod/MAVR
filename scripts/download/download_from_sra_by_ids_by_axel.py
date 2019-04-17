#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Collections.General import IdList


#from RouToolPa.Tools.Abstract import Tool
#from RouToolPa.Routines import NCBIRoutines

from RouToolPa.Tools.LinuxTools import Axel

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--ids", action="store", dest="ids",
                    type=lambda s: s.split(","),
                    help="Comma-separated list of SRA ids to download")
parser.add_argument("-f", "--id_file", action="store", dest="id_file",
                    help="File with SRA ids(one per line) to download")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of simultaneous downloads")
parser.add_argument("-c", "--connections", action="store", dest="connections", type=int, default=8,
                    help="Number of connections for each download")

args = parser.parse_args()

if (not args.ids) and (not args.id_file):
    raise ValueError("Both ids and id file were not set")

loader = IdList()
id_list = loader.read(args.id_file) if args.id_file else args.ids

Axel.threads = args.threads
Axel.parallel_download_from_sra(id_list, args.connections)
"""
options_list = []
for entry_id in id_list:
    ftp_path = NCBIRoutines.get_sra_ftp_path_from_id(entry_id)
    options_list.append("-n %i %s" % (args.connections, ftp_path))

tool = Tool(cmd="axel", max_threads=args.threads)

tool.parallel_execute(options_list)

for filename in os.listdir(os.getcwd()):
    if ".sra" not in filename:
        continue
    tool.safe_mkdir(filename[:-4])
    os.system("mv %s %s/" % (filename, filename[:-4]))
"""