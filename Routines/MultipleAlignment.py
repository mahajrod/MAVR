__author__ = 'mahajrod'
import os

from collections import OrderedDict
import numpy as np

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.Data import CodonTable
from Bio.Align import MultipleSeqAlignment

from Bio.SeqRecord import SeqRecord

from CustomCollections.GeneralCollections import SynDict, TwoLvlDict
from Routines.Sequence import SequenceRoutines


class MultipleAlignmentRoutines(SequenceRoutines):
    def __init__(self):
        SequenceRoutines.__init__(self)
        self.tmp_count_dict = {"A": 0, "T": 0, "G": 0, "C": 0, "N": 0}
        self.strict_nucleotide_list = ["A", "T", "G", "C"]
        self.strict_nucleotide_set = {"A", "T", "G", "C", "N"}
        self.nucleotide_list = ["A", "T", "G", "C", "N"]
        self.nucleotide_set = {"A", "T", "G", "C", "N"}

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
    def get_position_presence_matrix(alignment, gap_symbol="-", verbose=True):

        number_of_sequences = len(alignment)
        alignment_length = len(alignment[0])
        # converting alignment to numpy letter array stored by columns!
        alignment_array = np.array([list(rec) for rec in alignment], np.character, order="F")

        if verbose:
            print("%i sequences in alignment" % number_of_sequences)
            print("%i columns in alignment" % alignment_length)

        position_presence_array = np.array([[0 for letter in rec.seq] for rec in alignment], int, order="F")

        for column in range(0, alignment_length):
            column_list = alignment_array[:, column]
            number_of_gaps = 0
            for element in column_list:
                if element == gap_symbol:
                    number_of_gaps += 1
            for row in range(0, number_of_sequences):
                position_presence_array[row, column] = -number_of_gaps if alignment_array[row, column] == gap_symbol else number_of_sequences - number_of_gaps

        #for row in range(0, number_of_sequences):
        #    for column in range(0, alignment_length):
        #print alignment_array[0, ]
        #print alignment_array[:, 1]
        #print alignment_array
        #print position_presence_array
        return position_presence_array

    def get_position_presence_matrix_fom_file(self, alignment_file, output_file, format="fasta", gap_symbol="-",
                                              verbose=True):

        alignment = AlignIO.read(alignment_file, format=format)

        position_matrix = self.get_position_presence_matrix(alignment, gap_symbol, verbose=verbose)
        np.savetxt(output_file, position_matrix, fmt="%i", delimiter='\t')
        #print position_matrix
        return position_matrix

    def count_unique_positions_per_sequence(self, alignment, gap_symbol="-", verbose=True, ):
        number_of_sequences = len(alignment)
        alignment_length = len(alignment[0])
        position_presence_matrix = self.get_position_presence_matrix(alignment, gap_symbol=gap_symbol, verbose=verbose)
        unique_position_count_dict = OrderedDict()

        for row in range(0, number_of_sequences):
            sequence_id = alignment[row].id
            unique_positions = 0
            for column in range(0, alignment_length):
                if (position_presence_matrix[row, column] == 1) or (position_presence_matrix[row, column] == -1):
                    unique_positions += 1

            unique_position_count_dict[sequence_id] = unique_positions

        return unique_position_count_dict

    def count_unique_positions_per_sequence_from_file(self, alignment_file, output_prefix, format="fasta",
                                                      gap_symbol="-", return_mode="absolute", verbose=True):

        alignment = AlignIO.read(alignment_file, format=format)
        number_of_sequences = len(alignment)
        alignment_length = len(alignment[0])
        position_presence_matrix = self.get_position_presence_matrix(alignment, gap_symbol=gap_symbol, verbose=verbose)
        unique_position_count_dict = SynDict()
        unique_position_count_percent_dict = SynDict()

        for row in range(0, number_of_sequences):
            sequence_id = alignment[row].id
            unique_positions = 0
            for column in range(0, alignment_length):
                if (position_presence_matrix[row, column] == 1) or (position_presence_matrix[row, column] == -1):
                    unique_positions += 1

            unique_position_count_dict[sequence_id] = unique_positions
            unique_position_count_percent_dict[sequence_id] = 100 * float(unique_positions) / (alignment_length - str(alignment[row].seq).count(gap_symbol))

        unique_position_count_dict.write("%s.absolute_counts" % output_prefix)
        unique_position_count_percent_dict.write("%s.percent_counts" % output_prefix)

        return unique_position_count_dict if return_mode == "absolute" else unique_position_count_percent_dict

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
                print("Alignment was not filtered: allowed number of sequence with gaps is equal to number of sequences")
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

    def extract_degenerate_sites_from_codon_alignment(self, alignment, genetic_code_table=1,
                                                      remove_codon_columns_with_Ns=False):
        degenerate_codon_set = self.get_degenerate_codon_set(genetic_code_table)
        number_of_alignments = len(alignment)
        alignment_length = len(alignment[0])
        if alignment_length % 3 > 0:
            raise(ValueError, "Length of alignment is not divisible by 3")
        else:
            number_of_codons = int(alignment_length / 3)
        degenerate_columns = []
        degenerate_codons = []
        for i in range(0, number_of_codons):
            position_strings = []
            for j in range(0, 3):
                position_strings.append(list(set(alignment[:, 3*i + j])))
            if remove_codon_columns_with_Ns:
                if ("N" in position_strings[0]) or ("N" in position_strings[1]) or ("N" in position_strings[2]):
                    continue
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
                for k in 0, 1, 2:
                    degenerate_codons.append(alignment[:, 3*i + k])
        number_of_degenerate_columns = len(degenerate_columns)
        record_list = []
        degenerate_codons_record_list = []
        for i in range(0, number_of_alignments):
            string = ""
            for j in range(0, number_of_degenerate_columns):
                string += degenerate_columns[j][i]
            record = SeqRecord(seq=Seq(string), id=alignment[i].id, description=alignment[i].description)
            record_list.append(record)
            string = ""
            for j in range(0, 3 *number_of_degenerate_columns):
                string += degenerate_codons[j][i]
            record = SeqRecord(seq=Seq(string), id=alignment[i].id, description=alignment[i].description)
            degenerate_codons_record_list.append(record)


        #print(number_of_alignments, alignment_length)
        degenerate_alignment = MultipleSeqAlignment(record_list)
        degenerate_codons_alignment = MultipleSeqAlignment(degenerate_codons_record_list)
        return degenerate_alignment, degenerate_codons_alignment

    def extract_sites_by_column_expression_from_alignment(self, alignment, expression, remove_columns_with_Ns=False,
                                                          remove_columns_with_ambigous_nucleotides=False):
        number_of_alignments = len(alignment)
        alignment_length = len(alignment[0])
        extracted_columns = []

        for i in range(0, alignment_length):
            column_string = alignment[:, i]
            if remove_columns_with_Ns and ("N" in column_string):
                continue

            if remove_columns_with_ambigous_nucleotides:
                nucleotide_set = {s for s in column_string}
                if nucleotide_set - self.nucleotide_set:
                    continue
            if not expression(column_string):
                continue

            extracted_columns.append(column_string)

        number_of_extracted_columns = len(extracted_columns)
        record_list = []

        for i in range(0, number_of_alignments):
            string = ""
            for j in range(0, number_of_extracted_columns):
                string += extracted_columns[j][i]
            record = SeqRecord(seq=Seq(string), id=alignment[i].id, description=alignment[i].description)
            record_list.append(record)

        #print(number_of_alignments, alignment_length)
        extracted_alignment = MultipleSeqAlignment(record_list)

        return extracted_alignment

    def variable_alignment_column(self, column_string):
        different_nucleotide_number = 0

        for nucleotide in self.strict_nucleotide_list:
            self.tmp_count_dict[nucleotide] = column_string.count(nucleotide)
            if self.tmp_count_dict[nucleotide] > 0:
                different_nucleotide_number += 1

        return True if different_nucleotide_number > 1 else False

    def parsimony_informative_column(self, column_string):
        nucleotides_with_at_least_two_counts = 0

        for nucleotide in self.strict_nucleotide_list:
            self.tmp_count_dict[nucleotide] = column_string.count(nucleotide)
            if self.tmp_count_dict[nucleotide] >= 2:
                nucleotides_with_at_least_two_counts += 1

        return True if nucleotides_with_at_least_two_counts >= 2 else False

    def extract_variable_sites_from_alignment(self,
                                              alignment,
                                              remove_columns_with_Ns=False,
                                              remove_columns_with_ambigous_nucleotides=False):

        return self.extract_sites_by_column_expression_from_alignment(alignment,
                                                                      self.variable_alignment_column,
                                                                      remove_columns_with_Ns=remove_columns_with_Ns,
                                                                      remove_columns_with_ambigous_nucleotides=remove_columns_with_ambigous_nucleotides)

    def extract_variable_sites_from_alignment_from_file(self,
                                                        alignment_file,
                                                        output_file,
                                                        format="fasta",
                                                        remove_columns_with_Ns=False,
                                                        remove_columns_with_ambigous_nucleotides=False):

        alignment = AlignIO.read(alignment_file, format=format)
        variable_sites_alignment = self.extract_variable_sites_from_alignment(alignment,
                                                                              remove_columns_with_Ns=remove_columns_with_Ns,
                                                                              remove_columns_with_ambigous_nucleotides=remove_columns_with_ambigous_nucleotides)

        AlignIO.write([variable_sites_alignment], output_file, format=format)

    def extract_parsimony_informative_sites_from_alignment(self,
                                                           alignment,
                                                           remove_columns_with_Ns=False,
                                                           remove_columns_with_ambigous_nucleotides=False):

        return self.extract_sites_by_column_expression_from_alignment(alignment,
                                                                      self.parsimony_informative_column,
                                                                      remove_columns_with_Ns=remove_columns_with_Ns,
                                                                      remove_columns_with_ambigous_nucleotides=remove_columns_with_ambigous_nucleotides)

    def extract_parsimony_informative_sites_from_alignment_from_file(self,
                                                                     alignment_file,
                                                                     output_file,
                                                                     format="fasta",
                                                                     remove_columns_with_Ns=False,
                                                                     remove_columns_with_ambigous_nucleotides=False):

        alignment = AlignIO.read(alignment_file, format=format)
        variable_sites_alignment = self.extract_parsimony_informative_sites_from_alignment(alignment,
                                                                                           remove_columns_with_Ns=remove_columns_with_Ns,
                                                                                           remove_columns_with_ambigous_nucleotides=remove_columns_with_ambigous_nucleotides)

        AlignIO.write([variable_sites_alignment], output_file, format=format)

    def extract_degenerate_sites_from_codon_alignment_from_file(self, alignment_file, output_prefix,
                                                                genetic_code_table=1, format="fasta",
                                                                remove_codon_columns_with_Ns=False):
        alignment = AlignIO.read(alignment_file, format=format)
        degenerated_alignment, degenerate_codons_alignment = self.extract_degenerate_sites_from_codon_alignment(alignment,
                                                                                                                genetic_code_table=genetic_code_table,
                                                                                                                remove_codon_columns_with_Ns=remove_codon_columns_with_Ns)

        output_alignment_file = "%s.4fold_degenerated_sites.fasta" % output_prefix
        output_codon_alignment_file = "%s.codons_with_4fold_degenerated_sites.fasta" % output_prefix
        AlignIO.write([degenerated_alignment], output_alignment_file, format=format)
        AlignIO.write([degenerate_codons_alignment], output_codon_alignment_file, format=format)


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

    @staticmethod
    def translate_codon_alignment(codon_alignment_file, protein_alignment_file, format="fasta",
                                  gap_symbol="-", table=1):
        alignment = AlignIO.read(codon_alignment_file, format=format)

        def translated_record_generator(alignment):
            for record in alignment:
                yield SeqRecord(seq=record.seq.translate(gap=gap_symbol, table=table), id=record.id,
                                description=record.description)

        SeqIO.write(translated_record_generator(alignment), protein_alignment_file, format=format)

    def gather_alignment_stats(self, alignment, gap_symbol="-", verbose=False):

        number_of_sequences, alignment_length = self.get_general_statistics(alignment, verbose=verbose)

        unique_position_count_dict = self.count_unique_positions_per_sequence(alignment, gap_symbol=gap_symbol,
                                                                              verbose=verbose)
    @staticmethod
    def get_seq_id_dict(alignment):
        seq_id_index_dict = SynDict()
        for i in range(0, len(alignment)):
            seq_id_index_dict[i] = alignment[i].id
        return seq_id_index_dict

    def count_dNdS_by_reference_seq_in_codon_alignment(self, codon_alignment, reference_seq_id, genetic_code_table=1,
                                                       gap_symbol_list=["-",], use_ambigious_table=False,
                                                       output_file=None):

        number_of_sequences = len(codon_alignment)
        alignment_length = len(codon_alignment[0])
        if alignment_length % 3 > 0:
            raise ValueError("Length of alignment is not divisible by 3")
        else:
            number_of_codons = int(alignment_length / 3)

        seq_id_index_dict = self.get_seq_id_dict(codon_alignment)
        for i in range(0, number_of_sequences):
            if seq_id_index_dict[i] == reference_seq_id:
                reference_seq_index = i
                break
        else:
            raise ValueError("Reference sequence id (%s) is absent in alignment" % reference_seq_id)

        translation_table = CodonTable.ambiguous_dna_by_id[genetic_code_table] if use_ambigious_table else CodonTable.unambiguous_dna_by_id[genetic_code_table]
        dN_dS_W_dict = TwoLvlDict()

        for i in range(0, number_of_sequences):
            if i != reference_seq_index:
                dN_dS_W_dict[seq_id_index_dict[i]] = OrderedDict({"dN": 0, "dS": 0, "W": 'NA'})

        for i in range(0, number_of_codons):
            reference_codon = codon_alignment[reference_seq_index, i*3: (i+1) * 3]
            for symbol in gap_symbol_list:
                if symbol in reference_codon:
                    break
            else:
                reference_aminoacid = translation_table.forward_table[reference_codon]
                if reference_aminoacid == "X":
                    continue
                for j in range(0, number_of_sequences):
                    if j == reference_seq_index:
                        continue
                    sequence_codon = codon_alignment[j, i * 3: (i + 1) * 3]
                    for symbol in gap_symbol_list:
                        if symbol in sequence_codon:
                            break
                    else:
                        sequence_aminoacid = translation_table.forward_table[sequence_codon]
                        if sequence_aminoacid == "X":
                            continue

                        if (sequence_aminoacid == reference_aminoacid) and (sequence_codon != reference_codon):
                            dN_dS_W_dict[seq_id_index_dict[j]]["dS"] += 1
                        elif sequence_aminoacid != reference_aminoacid:
                            dN_dS_W_dict[seq_id_index_dict[j]]["dN"] += 1

        for sequence_id in dN_dS_W_dict:
            if dN_dS_W_dict[sequence_id]["dS"] != 0:
                dN_dS_W_dict[sequence_id]["W"] = float(dN_dS_W_dict[sequence_id]["dN"]) / float(dN_dS_W_dict[sequence_id]["dS"])

        if output_file:
            dN_dS_W_dict.write(output_file, absent_symbol="NA")

        return dN_dS_W_dict

    def count_dNdS_by_reference_seq_in_codon_alignment_from_file(self, alignment_file, reference_seq_id,
                                                                 genetic_code_table=1, gap_symbol_list=["-",],
                                                                 use_ambigious_table=False, output_file=None,
                                                                 format="fasta"):

        codon_alignment = self.parse_alignment(alignment_file, filetype=format)

        return self.count_dNdS_by_reference_seq_in_codon_alignment(codon_alignment, reference_seq_id,
                                                                   genetic_code_table=genetic_code_table,
                                                                   gap_symbol_list=gap_symbol_list,
                                                                   use_ambigious_table=use_ambigious_table,
                                                                   output_file=output_file)

    @staticmethod
    def extract_codon_positions(alignment):

        return [alignment[:, i::3] for i in 0, 1, 2]

    def extract_codon_positions_from_file(self, alignment_file, output_prefix, format="fasta"):
        alignment = self.parse_alignment(alignment_file, filetype="fasta")

        codon_position_alignments = self.extract_codon_positions(alignment)

        for position in 0, 1, 2:
            output_file = "%s.pos_%i%s" % (output_prefix, position + 1, self.split_filename(alignment_file)[-1])
            AlignIO.write(codon_position_alignments[position], output_file, format=format)


    """
    @staticmethod
    def parse_alignment(input_file, filetype="fasta"):
        print("Parsing alignment file %s" % input_file)
        #print(input_file)
        alignment = AlignIO.read(input_file, filetype)
        #print(alignment)
        print("Totaly %i sequences in alignment" % len(alignment))
        return alignment
    """