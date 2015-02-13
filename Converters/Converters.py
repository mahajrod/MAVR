#!/usr/bin/env python

import collections
import os
import re
from Bio import SeqIO
from Bio import AlignIO
from Bio import Phylo
from Bio.Alphabet import IUPAC, Gapped
from BCBio import GFF
from Parsers.General import parse_metamiga_fasta


def fastq2fasta(input_file, output_file):
    input_fd = open(input_file, "r")
    output_fd = open(output_file, "w")

    for line in input_fd:
        if line[0] == "@":
            output_fd.write(">" + line[1:])
            i = 0
            for line in input_fd:
                if line[0] == "+":
                    break
                i += 1
                output_fd.write(line)
            i -= 1
            for line in input_fd:
                if i == 0:
                    break
                i -= 1

    input_fd.close()
    output_fd.close()


def fasta2phylip(input_file, output_file):
    alignment = AlignIO.read(input_file, "fasta")
    number_of_sequences = len(alignment)
    length_of_alignment = alignment.get_alignment_length()
    record_id_list = []
    for record in alignment:
        record_id_list.append(len(record.id))
    max_id_length = max(record_id_list)
    fd = open(output_file, "w")
    fd.write("%i %i\n" % (number_of_sequences, length_of_alignment))
    for record in alignment:
        spaces = " " * (1 + max_id_length - len(record.id))
        fd.write(record.id + spaces + str(record.seq) + "\n")
    fd.close()


def fasta2nexus_TreeSAAP(input_file, output_file):
    #converts fasta to nexus for TreeSAAP
    alignment = AlignIO.read(input_file, "fasta")
    #number_of_sequences = len(alignment)
    #length_of_alignment = alignment.get_alignment_length()
    record_id_list = []
    #for record in alignment:
    #    record_id_list.append(len(record.id))
    #max_id_length = max(record_id_list)
    with open(output_file, "w") as fd:
        fd.write("#NEXUS\nbegin data;\n\tformat noninterleave datatype: DNA\nmatrix\n")
        for record in alignment:
            fd.write(record.id + "\n")
            fd.write(str(record.seq) + "\n")
        fd.write(";\nEND;\n")


def gff22gff3(input_file, output_file, target_lines=100000):

    in_fd = open(input_file, "r")
    out_fd = open(output_file, "w")
    GFF.write(GFF.parse(in_fd, target_lines=target_lines), out_fd)
    in_fd.close()
    out_fd.close()


def gff32gtf(input_file, output_file):
    os.system("gffread %s -T -o %s" % (input_file, output_file))


def gtf2gff3(input_file, output_file):
    os.system("gffread %s -o %s" % (input_file, output_file))


def convert_metamiga_fasta(input_fasta, output_fasta):
    parsed_fasta = parse_metamiga_fasta(input_fasta)
    #seq_list =
    #print(seq_list[0])
    #print("\n")
    #print(seq_list[0].seq)
    SeqIO.write(list(parsed_fasta.values()), output_fasta, "fasta")


def convert_alignment(input_file, input_filetype, output_file, output_filetype, alphabet=Gapped(IUPAC.ambiguous_dna)):
    alignment = AlignIO.parse(input_file, input_filetype, alphabet=alphabet)
    AlignIO.write(alignment, output_file, output_filetype)


def convert_tree(input_file, input_filetype, output_file, output_filetype):
    #tree = Phylo.read(input_file, input_filetype)
    Phylo.convert(input_file, input_filetype, output_file, output_filetype)
    #Phylo.write(tree, output_file, output_filetype)


def convert_sequences(input_file, input_index, input_filetype, output_file, output_filetype):
    #check if input is list-like object(checking for list methods in object)
    if isinstance(input_file, collections.MutableSequence):
        input_data = input_file
    else:
        input_data = [input_file]

    record_dict = SeqIO.index_db(input_index, input_data, input_filetype)
    SeqIO.write(record_dict.values(), output_file, output_filetype)

def sam2gtf(input_file, output_file, start_position=1, source="custom", scaffold=19):
    #TODO: works only for transcript alignments!!!!!!!!!!!!
    isoform = ["a", "b", "c", "d"]
    shift = start_position - 1
    isoform_index = 0
    with open(input_file, "r") as in_fd:
        with open(output_file, "w") as out_fd:
            for line in in_fd:
                if line[0] == "@":
                    continue
                sam_line = line.strip().split("\t")
                id = sam_line[0]
                region = sam_line[2]
                start = int(sam_line[3])
                quality = sam_line[4]
                cigar_string = sam_line[5]
                splited = re.split("M|N", cigar_string)[:-1]
                print(splited)
                segments_lengthes = list(map(lambda x: int(x), splited))
                isoform_name = "nxf1-" + isoform[isoform_index]
                isoform_index += 1
                transcript_line = '19\t%s\ttranscript\t%i\t%i\t.\t+\t.\tgene_id "NXF1"; gene_version "1"; transcript_id "%s"; transcript_version "1"; gene_name "nxf-1"; gene_source "custom"; gene_biotype "protein_coding"; transcript_name "%s";\n' % (source, start + shift, shift + sum(segments_lengthes), isoform_name, isoform_name)
                out_fd.write(transcript_line)
                for i in range(0, len(segments_lengthes)):
                    exon_start = shift + start + sum(segments_lengthes[:i])
                    if i % 2 == 0:
                        exon_line = '19\t%s\texon\t%i\t%i\t.\t+\t.\tgene_id "NXF1"; gene_version "1"; transcript_id "%s"; transcript_version "1"; exon_number "%i"; gene_name "nxf-1"; gene_source "custom"; gene_biotype "protein_coding"; transcript_name "%s";\n' % (source, exon_start, exon_start + segments_lengthes[i] - 1, isoform_name, i / 2, isoform_name)
                        out_fd.write(exon_line)



    #skip_header