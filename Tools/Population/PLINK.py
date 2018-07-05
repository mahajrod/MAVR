
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Tools.Abstract import Tool


class PLINK(Tool):
    def __init__(self, path="", max_threads=4, max_memory="100G", max_per_thread_memory="5G"):
        Tool.__init__(self, "plink", path=path, max_threads=max_threads, max_memory=max_memory, max_per_thread_memory=max_per_thread_memory)

        # constants for PLINK bed file
        self.samples_per_byte = 4
        self.signature_byte_number = 3

        # constans for PLINK bim file
        self.allel_columns_in_bim_file = (4, 5)

    @staticmethod
    def get_samples_list_from_plink_fam_file(plink_fam_file, verbose=False, sample_id_column=1):
        samples_list = []
        with open(plink_fam_file, "r") as plink_fam_fd:
            for line in plink_fam_fd:
                samples_list.append(line.split("\t")[sample_id_column])

        if verbose:
            print("%i samples were found" % len(samples_list))

        return samples_list

    def extract_sequences_from_plink_binary_snp_data(self, plink_binary_files_prefix, output_file, verbose=False,
                                                     output_format="fasta"):

        plink_bed_file = "%s.bed" % plink_binary_files_prefix
        plink_bim_file = "%s.bim" % plink_binary_files_prefix
        plink_fam_file = "%s.fam" % plink_binary_files_prefix

        if verbose:
            print("Parsing .fam file...")
        samples_list = self.get_samples_list_from_plink_fam_file(plink_fam_file, verbose=verbose)
        sample_number = len(samples_list)

        bytes_per_snp = sample_number / self.samples_per_byte if sample_number % self.samples_per_byte == 0 else int(sample_number/self.samples_per_byte) + 1

        if verbose:
            print("Parsing .bim file...")
        allels_array = np.loadtxt(plink_bim_file, usecols=self.allel_columns_in_bim_file, dtype=np.str)
        snp_number = len(allels_array)

        if verbose:
            print("Reading .bed file...")
        with open(plink_bed_file, 'rb') as plink_bed_fd:
            plink_bed_fd.read(self.signature_byte_number)
            genotypes_array = np.reshape(np.fromfile(plink_bed_fd, dtype='uint8'), (-1, bytes_per_snp))

        if snp_number != len(genotypes_array):
            raise ValueError("ERROR!!! .bed and .bim files are incompatible by length or files were parsed with errors!")

        record_list = []

        for sample_index in range(0, sample_number):
            if verbose:
                print("Handling sample %i (%s)..." % (sample_index + 1, samples_list[sample_index]))
            sample_byte = int(sample_index / self.samples_per_byte)
            sample_offset = sample_index % self.samples_per_byte

            sample_genotype_list = []
            for pos_index in xrange(0, snp_number):
                pos_genotype = (genotypes_array[pos_index, sample_byte] >> (2 * sample_offset)) & 0b11
                #print pos_index
                #print allels_array[pos_index]
                if pos_genotype == 0b00:
                    #print allels_array[pos_genotype][0]
                    sample_genotype_list.append(allels_array[pos_index][0])

                elif pos_genotype == 0b01:
                    #print "N"
                    sample_genotype_list.append("N")
                elif pos_genotype == 0b10:
                    sample_genotype_list.append(self.ambiguous_nucleotides_string_reverse_dict[allels_array[pos_index][0] + allels_array[pos_index][1]
                                                                                          if allels_array[pos_index][0] < allels_array[pos_index][1]
                                                                                          else allels_array[pos_index][1] + allels_array[pos_index][0]])
                    #print self.ambiguous_nucleotides_string_reverse_dict[allels_array[pos_index][0] + allels_array[pos_index][1]
                    #                                                                      if allels_array[pos_index][0] < allels_array[pos_index][1]
                    #                                                                      else allels_array[pos_index][1] + allels_array[pos_index][0]]
                else:
                    sample_genotype_list.append(allels_array[pos_index][1])
                    #print allels_array[pos_index][1]

            record_list.append(SeqRecord(id=samples_list[sample_index],
                                         description="",
                                         seq=Seq("".join(sample_genotype_list))))
            #if verbose:
            #    print record_list[sample_index]

        SeqIO.write(record_list, output_file, format=output_format)
