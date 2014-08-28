#!/usr/bin/env python

from Bio import AlignIO, SeqIO


def matched_record_generator(record_dict, white_list):
    for record in white_list:
        yield record_dict[record]


def white_list_generator(white_list_file):
    with open(white_list_file, "r") as in_fd:
        for line in in_fd:
            tmp = line.strip()
            if tmp != "":
                yield tmp

alignment_file = "/home/mahajrod/genetics/MH_selection/ExaML_tree/merged_no_rRNA.fasta"
"""
output_files = ["/home/mahajrod/genetics/MH_selection/ExaML_tree/alignment_taxa1.fasta",
                "/home/mahajrod/genetics/MH_selection/ExaML_tree/alignment_taxa2.fasta",
                "/home/mahajrod/genetics/MH_selection/ExaML_tree/alignment_taxa3.fasta"]
white_list_files = ["/home/mahajrod/genetics/MH_selection/ExaML_tree/part1.taxa",
                    "/home/mahajrod/genetics/MH_selection/ExaML_tree/part2.taxa",
                    "/home/mahajrod/genetics/MH_selection/ExaML_tree/part3.taxa"]
"""
output_files = ["/home/mahajrod/genetics/MH_selection/ExaML_tree/alignment_taxa1_1.fasta",
                "/home/mahajrod/genetics/MH_selection/ExaML_tree/alignment_taxa1_2.fasta",
                "/home/mahajrod/genetics/MH_selection/ExaML_tree/alignment_taxa1_3.fasta"
                ]
white_list_files = ["/home/mahajrod/genetics/MH_selection/ExaML_tree/part1_1.taxa",
                    "/home/mahajrod/genetics/MH_selection/ExaML_tree/part1_2.taxa",
                    "/home/mahajrod/genetics/MH_selection/ExaML_tree/part1_3.taxa"
                    ]


alignment_records = SeqIO.to_dict(SeqIO.parse(alignment_file, "fasta"))
for white_list_file, output_file in zip(white_list_files, output_files):
    SeqIO.write(matched_record_generator(alignment_records, white_list_generator(white_list_file)), output_file, "fasta")


