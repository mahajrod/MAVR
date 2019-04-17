#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import collections
from Bio import SeqIO
from Bio import AlignIO
from RouToolPa.Parsers.General import parse_metamiga_fasta
from RouToolPa.Routines.Sequence import SequenceRoutines




class SequenceConverters(SequenceRoutines):

    def __init__(self):
        SequenceRoutines.__init__(self)

    @staticmethod
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

    @staticmethod
    def fasta2phylip(input_file, output_file, mode="naming_same"):
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
            if mode == "naming_same":
                spaces = " " * (2 + max_id_length - len(record.id))
                fd.write(record.id + spaces + str(record.seq) + "\n")
            else:
                fd.write(record.id + "\n")
                fd.write(str(record.seq) + "\n")
        fd.close()

    @staticmethod
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

    @staticmethod
    def convert_metamiga_fasta(input_fasta, output_fasta):
        parsed_fasta = parse_metamiga_fasta(input_fasta)
        #seq_list =
        #print(seq_list[0])
        #print("\n")
        #print(seq_list[0].seq)
        SeqIO.write(list(parsed_fasta.values()), output_fasta, "fasta")

    def convert_sequences(self, input_file, input_filetype,
                          output_file, output_filetype,
                          parsing_mode="parse", input_index="tmp.idx", alphabet=None):

        #check if input is list-like object(checking for list methods in object)
        if isinstance(input_file, collections.MutableSequence):
            input_data = input_file
        else:
            input_data = [input_file]

        alphab = eval("IUPAC.%s" % alphabet) if alphabet is not None else None

        #record_dict = self.parse_seq_file(input_data, parsing_mode, format=input_filetype, index_file=input_index)#SeqIO.index_db(input_index, input_data, input_filetype)
        #SeqIO.write(record_dict.values(), output_file, output_filetype, alphabet=alphab)
        SeqIO.convert(input_file, input_filetype, output_file, output_filetype, alphabet=alphab)
        if parsing_mode == "index_db":
            os.remove(input_index)

    def sequence2fastq(self, sequence_file, output_file, mode, default_quality,
                       format="fasta", index_file=None, score_type=33):
        #input_file, mode, format="fasta", index_file=None
        record_dict = self.parse_seq_file(sequence_file, mode, format=format, index_file=index_file)

        with open(output_file, "w") as out_fd:
            for record_id in record_dict:
                record_len = len(record_dict[record_id].seq)
                out_fd.write("@%s\n" % record_id)
                out_fd.write("%s\n+\n" % str(record_dict[record_id].seq))
                out_fd.write("%s\n" % (chr(score_type + default_quality) * record_len))
