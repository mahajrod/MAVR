__author__ = 'mahajrod'
import os

from Routines.SequenceCluster import SequenceClusterRoutines

from CustomCollections.GeneralCollections import SynDict, IdList


class GORoutines(SequenceClusterRoutines):
    def __init__(self):
        SequenceClusterRoutines.__init__(self)

    @staticmethod
    def extract_entries_by_GO_from_eggnogmapper_output(eggnogmapper_output, GO_file, output_prefix,
                                                       comments_prefix="#", separator="\t",
                                                       ):

        GO_list = IdList(filename=GO_file, column_number=0)

        #print "GOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO"
        #print GO_list
        print len(GO_list)

        extracted_entries_file = "%s.annotations" % output_prefix

        extracted_entries = 0

        with open(eggnogmapper_output, "r") as eggnog_fd:
            with open(extracted_entries_file, "w") as out_fd:
                for line in eggnog_fd:
                    if line[0] == comments_prefix:
                        out_fd.write(line)
                        continue
                    line_list = line.strip().split(separator)
                    entry_GO_list = line_list[5].split(",")
                    #print entry_GO_list
                    for GO in entry_GO_list:
                        if GO in GO_list:
                            out_fd.write(line)
                            extracted_entries += 1

        print("Extracted %i entries" % extracted_entries)