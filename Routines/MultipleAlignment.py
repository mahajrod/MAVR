__author__ = 'mahajrod'
import os

import numpy as np

from Bio.Seq import Seq
from Bio import SeqIO, AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

from Routines import SequenceRoutines
from CustomCollections.GeneralCollections import SynDict


class MultipleAlignmentRoutines:
    def __init__(self):
        pass

    @staticmethod
    def get_general_statistics(alignment, verbose=False):
        number_of_sequences = len(alignment)
        length_of_alignment = alignment.get_alignment_length()
        if verbose:
            print ("Alignment statistics:")
            print ("\tNumber of sequences %i" % number_of_sequences)
            print ("\tLength %i" % length_of_alignment)
        return number_of_sequences, length_of_alignment

    @staticmethod
    def get_position_presence_matrix(alignment, gap_symbol="-"):

        number_of_sequences = len(alignment)
        # converting alignment to numpy letter array stored by columns!
        align_array = np.array([list(rec) for rec in alignment], np.character, order="F")

        position_presence_array = np.array([0 for rec in alignment], int, order="F")
        print align_array[0, ]
        print align_array[:, 1]
        return align_array

    def get_position_presence_matrix_fom_file(self, alignment_file, output_file, format="fasta",
                                              gap_symbol="-"):
        alignment = AlignIO.read(alignment_file, format=format)

        position_matrix = self.get_position_presence_matrix(alignment, gap_symbol)

        print position_matrix

    @staticmethod
    def parse_alignment(input_file, filetype="fasta"):
        print("Parsing alignment file %s" % input_file)
        #print(input_file)
        alignment = AlignIO.read(input_file, filetype)
        #print(alignment)
        print("Totaly %i sequences in alignment" % len(alignment))
        return alignment

    @staticmethod
    def make_regions_from_columns_to_retain(columns_to_retain,):
        region_coordinates_list = [[columns_to_retain[0], columns_to_retain[0] + 1]]

        if len(columns_to_retain) == 1:
            return region_coordinates_list

        for i in range(1, len(columns_to_retain)):
            if region_coordinates_list[-1][1] == columns_to_retain[i]:
                region_coordinates_list[-1][1] += 1
            else:
                region_coordinates_list.append([columns_to_retain[i], columns_to_retain[i] + 1])

        return region_coordinates_list

    def remove_columns_with_gaps(self, alignment, maximum_number_of_gaps_in_column, gap_symbol="-", verbose=False):
        num_of_sequences, len_of_alignment = self.get_general_statistics(alignment, verbose=verbose)
        if maximum_number_of_gaps_in_column > num_of_sequences:
            raise ValueError("Allowed number of sequences with gap is bigger than number of sequences")
        if maximum_number_of_gaps_in_column == num_of_sequences:
            if verbose:
                print "Alignment was not filtered: allowed number of sequence with gaps is equal to number of sequences"
            return alignment

        columns_without_gaps_list = []
        for i in range(0, len_of_alignment):
            if alignment[:, i].count(gap_symbol) <= maximum_number_of_gaps_in_column:
                columns_without_gaps_list.append(i)

        region_to_retain = self.make_regions_from_columns_to_retain(columns_without_gaps_list)
        filtered_alignment = alignment[:, region_to_retain[0][0]: region_to_retain[0][1]]
        if len(region_to_retain) == 1:
            return filtered_alignment

        for start, end in region_to_retain[1:]:
            filtered_alignment += alignment[:, start: end]

        return filtered_alignment

    @staticmethod
    def slice_fasta_alignment(alignment_file, output_file, start, end):
        alignment = AlignIO.read(alignment_file, format="fasta")
        SeqIO.write(alignment[start:end], output_file, "fasta")

    @staticmethod
    def prepare_multigene_alignment_for_codeml(alignment_file, coordinates_file, output_file,
                                               format="fasta"):
        sequence_dict = SeqIO.to_dict(SeqIO.parse(alignment_file, format=format))
        number_of_sequences = len(sequence_dict)
        alignment_length = len(sequence_dict[sequence_dict.keys()[0]])

        gene_lengths_in_codons = np.loadtxt(coordinates_file, dtype=int, comments='#', usecols=(0,))/3
        number_of_genes = len(gene_lengths_in_codons)
        codon_number_string = " ".join(map(str, gene_lengths_in_codons))
        with open(output_file, "w") as out_fd:
            out_fd.write("%i %i G\n" % (number_of_sequences, alignment_length))
            out_fd.write("G %i %s\n" % (number_of_genes, codon_number_string))
            for record_id in sequence_dict:
                out_fd.write("%s\n%s\n" % (record_id, str(sequence_dict[record_id].seq)))

    @staticmethod
    def get_db_ids(search_dict):
        id_set = set()
        for query_id in search_dict:
            for hit in search_dict[query_id]:
                id_set.add(hit.id)
        return id_set

    @staticmethod
    def get_codon_alignment(protein_alignment, nucleotide_seq_dict, codon_alignment_file,
                            protein2cds_accordance_dict=None):
        codon_alignment = {}
        if protein2cds_accordance_dict:
            for record in protein_alignment:
                if record.id not in protein2cds_accordance_dict:
                    print("Corresponding CDS was not found for %s protein" % record.id)
                    return -1
        for record in protein_alignment:
            nucleotide_seq = ""
            i = 0
            for aminoacid in record.seq:
                if aminoacid == "-":
                    nucleotide_seq += "---"
                    continue
                else:
                    nucleotide_seq += str(nucleotide_seq_dict[protein2cds_accordance_dict[record.id] if protein2cds_accordance_dict else record.id].seq[3*i:3*(i+1)])
                    i += 1
            codon_alignment[record.id] = SeqRecord(Seq(nucleotide_seq),
                                                   id=record.id,
                                                   description=record.description,
                                                   name=record.name)
            #print(record.id, record.seq)
        SeqIO.write(list(codon_alignment.values()), codon_alignment_file, "fasta")
        return codon_alignment

    def get_codon_alignment_from_files(self, protein_aln_file, nucleotide_seq_file, codon_alignment_file,
                                       cds2protein_accordance_file=None,
                                       alignment_format="fasta", nucleotide_sequence_format="fasta",
                                       cds_index_file=None, retain_cds_index=False):
        protein_aln_dict = AlignIO.read(protein_aln_file, format=alignment_format)
        nucleotide_seq_dict = SeqIO.index_db(cds_index_file if cds_index_file else "nuc_tmp.idx", nucleotide_seq_file,
                                             format=nucleotide_sequence_format)

        protein2cds_accordance_dict = None
        if cds2protein_accordance_file:
            protein2cds_accordance_dict = SynDict()
            protein2cds_accordance_dict.read(cds2protein_accordance_file, key_index=1, value_index=0)

        self.get_codon_alignment(protein_aln_dict, nucleotide_seq_dict, codon_alignment_file,
                                 protein2cds_accordance_dict=protein2cds_accordance_dict)
        if (not cds_index_file) and (not retain_cds_index):
            os.remove("nuc_tmp.idx")

    @staticmethod
    def merge_alignment(alignment_file_list, merged_alignment_file, coordinates_file, format="fasta"):
        #print("Merging alignments...")
        alignment_list = []
        sequence_lengthes = []
        #print(alignment_file_list)

        alignment_file_list_sorted = sorted(alignment_file_list)
        #print(alignment_file_list_sorted)
        for alignment_file in alignment_file_list_sorted:
            #alignment_file.sort()
            parsed = AlignIO.read(alignment_file, format=format)
            parsed.sort()
            alignment_list.append(parsed)
        merged_alignment = None
        for alignment in alignment_list:
            if not merged_alignment:
                sequence_lengthes.append(alignment.get_alignment_length())
                merged_alignment = alignment
                continue
            #print(alignment)
            sequence_lengthes.append(alignment.get_alignment_length())
            merged_alignment += alignment
        SeqIO.write(merged_alignment, merged_alignment_file, "fasta")
        sequence_coordinates = []
        #
        for seq_length in sequence_lengthes:
            if not sequence_coordinates:
                sequence_coordinates.append((1, seq_length))
                continue
            sequence_coordinates.append((sequence_coordinates[-1][1]+1, sequence_coordinates[-1][1]+seq_length))
        #print(sequence_coordinates)
        with open(coordinates_file, "w") as coord_fd:
            coord_fd.write("#length\tstart\tend\n")
            for coord_tuple in sequence_coordinates:
                coord_fd.write("%i\t%i\t%i\n" % (coord_tuple[1] - coord_tuple[0] + 1, coord_tuple[0], coord_tuple[1]))
        return merged_alignment, sequence_lengthes, sequence_coordinates

    @staticmethod
    def extract_degenerate_sites_from_codon_alignment(alignment, genetic_code_table=1):
        degenerate_codon_set = SequenceRoutines.get_degenerate_codon_set(genetic_code_table)
        number_of_alignments = len(alignment)
        alignment_length = len(alignment[0])
        if alignment_length % 3 > 0:
            raise(ValueError, "Length of alignment is not divisible by 3")
        else:
            number_of_codons = int(alignment_length / 3)
        degenerate_columns = []
        for i in range(0, number_of_codons):
            position_strings = []
            for j in range(0, 3):
                position_strings.append(list(set(alignment[:, 3*i + j])))
            if (len(position_strings[0]) > 1) or (len(position_strings[1]) > 1):
                continue
            ambigious_codon = position_strings[0][0] + position_strings[1][0] + "N"
            """
            if Seq(ambigious_codon).translate(table=genetic_code_table) == "X":
                continue
            else:
                degenerate_columns.append(alignment[:, 3*i + 2])
                #print(i*3 +3)
            """
            if ambigious_codon in degenerate_codon_set:
                degenerate_columns.append(alignment[:, 3*i + 2])

        number_of_degenerate_columns = len(degenerate_columns)
        record_list = []
        for i in range(0, number_of_alignments):
            string = ""
            for j in range(0, number_of_degenerate_columns):
                string += degenerate_columns[j][i]
            record = SeqRecord(seq=Seq(string), id=alignment[i].id)
            record_list.append(record)

        #print(number_of_alignments, alignment_length)
        degenerate_alignment = MultipleSeqAlignment(record_list)

        return degenerate_alignment

    def extract_degenerate_sites_from_codon_alignment_from_file(self, alignment_file, output_alignment_file,
                                                                genetic_code_table=1, format="fasta"):
        alignment = AlignIO.read(alignment_file, format=format)
        degenerate_alignment = self.extract_degenerate_sites_from_codon_alignment(alignment,
                                                                                  genetic_code_table=genetic_code_table)

        AlignIO.write([degenerate_alignment], output_alignment_file, format=format)

    @staticmethod
    def sequences_from_alignment_generator(alignments, gap_symbol="-"):
        for record in alignments:
            record.seq = Seq(str(record.seq).replace(gap_symbol, ""))
            yield record

    def extract_sequences_from_alignment(self, alignment_file, output_file, alignment_format="fasta",
                                         output_format="fasta", gap_symbol="-"):
        alignments = AlignIO.read(alignment_file, format=alignment_format)
        SeqIO.write(self.sequences_from_alignment_generator(alignments, gap_symbol=gap_symbol),
                    output_file, format=output_format)


