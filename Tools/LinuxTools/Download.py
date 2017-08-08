#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Routines import NCBIRoutines

from Tools.Abstract import Tool


class Wget(Tool):
    def __init__(self, path="", max_threads=4):

        Tool.__init__(self, "wget", path=path, max_threads=max_threads)

    def parallel_download_from_sra(self, id_list):
        options_list = []
        for entry_id in id_list:
            ftp_path = NCBIRoutines.get_sra_ftp_path_from_id(entry_id)
            options_list.append("--no-host-directories -rc -t 500 %s" % ftp_path)

        self.parallel_execute(options_list)

        for filename in os.listdir(os.getcwd()):
            if ".sra" not in filename:
                continue
            self.safe_mkdir(filename[:-4])
            os.system("mv %s %s/" % (filename, filename[:-4]))


class Axel(Tool):
    def __init__(self, path="", max_threads=4):

        Tool.__init__(self, "axel", path=path, max_threads=max_threads)

    def parallel_download_from_sra(self, id_list, connections=4):

        options_list = []
        for entry_id in id_list:
            ftp_path = NCBIRoutines.get_sra_ftp_path_from_id(entry_id)
            options_list.append("-n %i %s" % (connections, ftp_path))

        self.parallel_execute(options_list)

        for filename in os.listdir(os.getcwd()):
            if ".sra" not in filename:
                continue
            self.safe_mkdir(filename[:-4])
            os.system("mv %s %s/" % (filename, filename[:-4]))

if __name__ == "__main__":
    pass
