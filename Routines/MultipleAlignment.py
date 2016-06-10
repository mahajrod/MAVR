__author__ = 'mahajrod'

from Bio import SeqIO
from Bio import AlignIO

from numpy import loadtxt

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

        gene_lengths_in_codons = loadtxt(coordinates_file, dtype=int, comments='#', usecols=(0,))/3
        number_of_genes = len(gene_lengths_in_codons)
        codon_number_string = " ".join(map(str, gene_lengths_in_codons))
        with open(output_file, "w") as out_fd:
            out_fd.write("%i %i G\n" % (number_of_sequences, alignment_length))
            out_fd.write("G %i %s\n" % (number_of_genes, codon_number_string))
            for record_id in sequence_dict:
                out_fd.write("%s\n%s\n" % (record_id, str(sequence_dict[record_id].seq)))




