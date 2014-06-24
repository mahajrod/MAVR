#!/usr/bin/env python
from Bio import SeqIO
from Bio import AlignIO
from Parse.ParseGeneral import parse_metamiga_fasta


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


def convert_metamiga_fasta(input_fasta, output_fasta):
    parsed_fasta = parse_metamiga_fasta(input_fasta)
    #seq_list =
    #print(seq_list[0])
    #print("\n")
    #print(seq_list[0].seq)
    SeqIO.write(list(parsed_fasta.values()), output_fasta, "fasta")