#!/usr/bin/env python
import os

from Tools.Abstract import Tool
from CustomCollections.GeneralCollections import IdSet


class RepeatMasker(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "repeatmasker", path=path, max_threads=max_threads)

    @staticmethod
    def convert_rm_out_to_gff(input_file, output_file, annotated_repeat_classes_file, annotated_repeat_families_file):
        repeat_classes_set = IdSet()
        repeat_families_set = IdSet()
        with open(input_file, "r") as in_fd:
            for i in range(0, 3):
                in_fd.readline()

            with open(output_file, "w") as out_fd:
                for line in in_fd:
                    tmp = line.strip().split()
                    strand = "+" if tmp[8] == "+" else "-"
                    repeat_class_family = tmp[10].split("/")
                    if len(repeat_class_family) == 1:
                        repeat_class_family.append(".")
                    repeat_classes_set.add(repeat_class_family[0])
                    repeat_families_set.add(repeat_class_family)
                    parameters = "Class=%s;Family=%s;Matching_repeat=%s;SW_score=%s;Perc_div=%s;Perc_del=%s;Pers_ins=%s" \
                                 % (repeat_class_family[0], repeat_class_family[1],
                                    tmp[9], tmp[0], tmp[1], tmp[2], tmp[3])
                    out_fd.write("%s\tRepeatMasker\trepeat\t%s\t%s\t.\t%s\t.\t%s\n"
                                 % (tmp[4], tmp[5], tmp[6], strand, parameters))
        repeat_classes_set.write(annotated_repeat_classes_file)
        repeat_families_set.write(annotated_repeat_families_file)
        #if annotated_repeat_types_file:
        """
        sed_string = "sed -r 's/.*Class=(.*);Family.*/\1/' %s | sort | uniq > %s" % (output_file,
                                                                                     annotated_repeat_classes_file)
        os.system(sed_string)
        """
    @staticmethod
    def extract_annotated_repeat_types_from_gff(gff_file, annotated_repeat_classes_file):
        sed_string = "sed -r 's/.*Class=(.*);Family.*/\1/' %s | sort | uniq > %s" % (gff_file,
                                                                                     annotated_repeat_classes_file)
        os.system(sed_string)

