#!/usr/bin/env python
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def parse_metamiga_fasta(input_file):
    print("Parsing %s" % input_file)
    seq_dict = {}
    fd = open(input_file, "r")
    sequence = ""
    for line in fd:
        if line[0] == ">":
            if sequence != "":
                seq_dict[seq_id] = SeqRecord(Seq(sequence), id=seq_id, description=gene, name=gene)
                seq_dict[seq_id].annotations["organism"] = species
            annotations_list = line[1:].strip().split("-")
            seq_id = annotations_list[0]
            species = annotations_list[1]
            gene = annotations_list[2]
            sequence = ""
        else:
            sequence += line.strip()
    seq_dict[seq_id] = SeqRecord(Seq(sequence), id=seq_id, description=gene, name=gene)
    seq_dict[seq_id].annotations["organism"] = species
    #print()
    fd.close()
    print("Totaly %i sequences were parsed" % len(seq_dict))
    #for seq in seq_dict:
    #    print(seq, seq_dict[seq])
    return seq_dict


def parse_alignment(input_file, filetype="fasta"):
    print("Parsing alignment file %s" % input_file)
    #print(input_file)
    alignment = AlignIO.read(input_file, filetype)
    #print(alignment)
    print("Totaly %i sequences in alignment" % len(alignment))
    return alignment


def parse_sv(input_file, separator="\t"):
    record_list = []
    with open(input_file, "r") as sv_fd:
        header = sv_fd.readline().strip()
        if header[0] == "#":
            header = header[1:]
        header = header.split("\t")
        for line in sv_fd:
            if line == "\n":
                break
            record_list.append(line.strip().split(separator))
    return header, record_list

"""
def parse_gene_list_file(input_file):
    #TODO: rewrite as class ?
    header, record_list = parse_sv()
"""