#!/usr/bin/env python2

#import sys
#import re
import os
#import glob
import argparse
from copy import deepcopy


class ScriptArg(object):
    """docstring for ScriptArg
        Parser cmd line arguments using argparse"""
    def __init__(self):
        #super(ScriptArg, self).__init__()
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument('-s', action='store', dest='source', help='Source files. If using more then one file, they should be comma separated')
        self.parser.add_argument('-g', action='store_true', dest='gui', help='Enables GUI')
        self.parser.add_argument('-v', action='store_true', dest='verbose', help='Enables extended output')
        self.arg = self.parser.parse_args()
        if self.arg.verbose:
            print("\n\tArguments of the script:\n\t\tSource: \t%s \n\t\tGUI: \t\t%s \n\t\tVerbose: \t%s\n" %
                  (self.arg.source, self.arg.gui, self.arg.verbose))


class SourceFiles(object):
    """docstring for SourceFiles"""
    ignore_list_files_ext = {
                            "index": [".idx"]
                            }
    file_types_ext = {
                    "genbank":[".genbank", ".gb"],
                    "fasta":[".fa", ".fasta"],
                    "vcf": [".vcf"],
                    "fastq": [".fastq", ".fq"],
                    "gff": [".gff"],
                    "gtf": [".gtf"]
                    }

    files_dict = {}
    ignored_files = []
    unrecognized_files = []
    filelist = []

    def __recognize_filetype(self, filename):
        file_ext = "." + filename.split(".")[-1]
        for known_filetype in self.file_types_ext:
            if file_ext in self.file_types_ext[known_filetype]:
                return known_filetype
        return None

    def filter_files_by_ext(self, black_list=[], white_list=[]):
        i = 0
        while i < len(self.filelist):
            file_type = self.__recognize_filetype(self.filelist[i])
            if not file_type:
                self.unrecognized_files.append(deepcopy(self.filelist[i]))
                self.filelist.pop(i)
                continue
            file_ext = "." + self.filelist[i].split(".")[-1]
            if black_list:
                if file_ext in black_list:
                    self.ignored_files.append(deepcopy(self.filelist[i]))
                    self.filelist.pop(i)
                    continue
            elif white_list:
                if file_ext not in black_list:
                    self.ignored_files.append(deepcopy(self.filelist[i]))
                    self.filelist.pop(i)
                    continue
            if file_type not in self.files_dict:
                self.files_dict[file_type] = [self.filelist[i]]
            else:
                self.files_dict[file_type].append(self.filelist[i])
            i += 1

        assert self.filelist, \
            "\tError 4: No files left after filtering"

    def __init__(self, source, verbose=False):
        #super(SourceFiles, self).__init__()
        self.source = source.split(",")
        self.verbose = verbose
        for source_path in self.source:
            assert os.path.exists(source_path), \
                "\tError 1: Source file/directory doesnt exist"
        for entry in self.source:
            if os.path.isfile(entry):
                self.filelist.append(entry)
            else:
                if entry[-1] != "/":			# to be sure that dirname ends with "/"
                    directory = entry + "/"
                else:
                    directory = entry

                dir_list = os.listdir(directory)
                #transform short filename to full
                for file_from_dir in dir_list:
                    self.filelist.append(directory + file_from_dir)
                #remove ignored files from list
        self.filter_files_by_ext()

        if self.verbose:
            print ("\tSource files:")
            for entry in self.filelist:
                print ("\t\t" + entry)
            #print ()
            if self.ignored_files:
                print ("\tIgnored files:")
                for entry in self.ignored_files:
                    print ("\t\t" + entry)
                #print ()
    """
    def sort_files(self):
        self.unrecognized_files = list(self.filelist)
        for entry in self.filelist:
            for ftype in self.file_types_ext:
                for extention in self.file_types_ext[ftype]:
                    #if entry[-self.max_ext_length:].find(extention) > -1 : #check for file type extension in last <max_ext_length> symbols of filename
                    if self.__get_file_ext(entry) == extention : #check for file type extension in last <max_ext_length> symbols of filename
                        if ftype not in self.files_dict: # check for presense filetype in output files_dict
                            self.files_dict[ftype] = []
                        if entry not in self.files_dict[ftype]:    # check for duplicates
                            self.unrecognized_files.remove(entry)
                            self.files_dict[ftype].append(entry)

        assert self.files_dict, \
            "\tError 3: None of files was recognized"

        if self.unrecognized_files:
            print ("\tWARNING: Next file(s) was(were) not recognized:")
            for entry in self.unrecognized_files:
                print ("\t\t" + entry)
        if self.verbose:
            print ("\tRecognized files:")
            for ftype in self.files_dict:
                print ("\t\t" + ftype + " format:")
                for files in self.files_dict[ftype]:
                    print ("\t\t\t" + files)
            print ()
        return self.files_dict
    """


def check_path(path_to_check):
    #returns path with / at end or blank path
    if path_to_check != "":
        if path_to_check[-1] != "/":
            return path_to_check + "/"
    return path_to_check
