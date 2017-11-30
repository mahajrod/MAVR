#!/usr/bin/env python
from Tools.Abstract import Tool
from CustomCollections.GeneralCollections import SynDict


class CDHit(Tool):
    """
    T
    """

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cd-hit", path=path, max_threads=max_threads)

    @staticmethod
    def convert_clustering_to_fam(input_clustering_file, output_fam):
        clusering_dict = SynDict()
        with open(input_clustering_file, "r") as in_fd:
            cluster_id = in_fd.readline()[1:].strip().replace(" ", "_")
            seq_list = []
            for line in in_fd:
                if line[0] == ">":
                    clusering_dict[cluster_id] = seq_list
                    cluster_id = in_fd.next()[1:].strip().replace(" ", "_")
                    seq_list = []
                    continue
                seq_list.append(line.strip().split()[2][1:-3])

        clusering_dict.write(output_fam, splited_values=True)



