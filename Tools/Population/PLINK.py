
import numpy as np
from collections import OrderedDict

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from Tools.Abstract import Tool
from Parsers.PLINK import PLINKReport

class PLINK(Tool):
    def __init__(self, path="", max_threads=4, max_memory="100G", max_per_thread_memory="5G"):
        Tool.__init__(self, "plink", path=path, max_threads=max_threads, max_memory=max_memory, max_per_thread_memory=max_per_thread_memory)

        # constants for PLINK bed file
        self.samples_per_byte = 4
        self.signature_byte_number = 3

        # constans for PLINK bim file
        self.allel_columns_in_bim_file = (4, 5)


    def parse_common_options(self, output_prefix, input_vcf_file=None, allow_noncanonical_chromosome_names=None, keep_autoconverted_files=None):

        options = " --threads %i" % self.threads
        options += " --out %s" % output_prefix
        options += " --allow-extra-chr" if allow_noncanonical_chromosome_names else ""
        options += " --keep-autoconv" if keep_autoconverted_files else ""
        options += " --vcf %s" % input_vcf_file if input_vcf_file else ""

        return options

    def find_runs_of_homozygosity_options(self,
                                          roh_calling_method=None,
                                          window_length_in_kb=None,
                                          min_homozygous_snps_per_window=None,
                                          max_heterozygous_snps_per_window=None,
                                          max_missing_snps_per_window=None,
                                          max_inverse_density_of_homozygous_snps_in_kb_per_snp=None,
                                          max_internal_gap_in_kb=None,
                                          min_roh_length=None,
                                          min_homozygous_snps_in_roh=None,
                                          min_scanning_window_hit_rate=None,
                                          generate_overlapping_segments=False,
                                          max_heterozygous_snps=None,
                                          min_concordance_across_jointly_homozygous_variants=None,
                                          homozygous_verbose=False):
        """
        ROH calling methods:
        <group | group-verbose> <consensus-match> <extend> <subtract-1-from-lengths>

        """
        options = " --homozyg %s" % roh_calling_method if roh_calling_method else " --homozyg"
        options += " --homozyg-window-snp %i" % min_homozygous_snps_per_window if min_homozygous_snps_per_window else ""                      # 50
        options += " --homozyg-window-het %i" % max_heterozygous_snps_per_window if max_heterozygous_snps_per_window else ""                  # 1
        options += " --homozyg-window-missing %i" % max_missing_snps_per_window if max_missing_snps_per_window else ""                        # 5
        options += " --homozyg-window-kb %i" % window_length_in_kb if window_length_in_kb else ""                                             # 5000
        options += " --homozyg-window-threshold %f" % min_scanning_window_hit_rate if min_scanning_window_hit_rate else ""                    # 0.05

        options += " --homozyg-density %i" % max_inverse_density_of_homozygous_snps_in_kb_per_snp if max_inverse_density_of_homozygous_snps_in_kb_per_snp else "" # 50
        options += " --homozyg-gap %i" % max_internal_gap_in_kb if max_internal_gap_in_kb else ""                                                                 # 1000
        options += " --homozyg-snp %i" % min_homozygous_snps_in_roh if min_homozygous_snps_in_roh else ""                                                         # 100
        options += " --homozyg-kb %i" % min_roh_length if min_roh_length else ""                                                                                  # 1000
        options += " --homozyg-match %f" % min_concordance_across_jointly_homozygous_variants if min_concordance_across_jointly_homozygous_variants else ""       # 0.95
        options += " --homozyg-het %i" % max_heterozygous_snps if max_heterozygous_snps else ""                                                                 # unlimited
        options += " --homozyg-group" if generate_overlapping_segments else ""
        options += " --homozyg-verbose" if homozygous_verbose else "" #

        return options

    def find_runs_of_homozygosity(self,
                                  output_prefix,
                                  input_vcf_file=None,
                                  allow_noncanonical_chromosome_names=None,
                                  keep_autoconverted_files=None,
                                  roh_calling_method=None,
                                  window_length_in_kb=None,
                                  min_homozygous_snps_per_window=None,
                                  max_heterozygous_snps_per_window=None,
                                  max_missing_snps_per_window=None,
                                  max_inverse_density_of_homozygous_snps_in_kb_per_snp=None,
                                  max_internal_gap_in_kb=None,
                                  min_roh_length=None,
                                  min_homozygous_snps_in_roh=None,
                                  min_scanning_window_hit_rate=None,
                                  generate_overlapping_segments=False,
                                  max_heterozygous_snps=None,
                                  min_concordance_across_jointly_homozygous_variants=None,
                                  homozygous_verbose=False):
        options = self.parse_common_options(output_prefix,
                                            input_vcf_file=input_vcf_file,
                                            allow_noncanonical_chromosome_names=allow_noncanonical_chromosome_names,
                                            keep_autoconverted_files=keep_autoconverted_files)

        options += self.find_runs_of_homozygosity_options(roh_calling_method=roh_calling_method,
                                                          window_length_in_kb=window_length_in_kb,
                                                          min_homozygous_snps_per_window=min_homozygous_snps_per_window,
                                                          max_heterozygous_snps_per_window=max_heterozygous_snps_per_window,
                                                          max_missing_snps_per_window=max_missing_snps_per_window,
                                                          max_inverse_density_of_homozygous_snps_in_kb_per_snp=max_inverse_density_of_homozygous_snps_in_kb_per_snp,
                                                          max_internal_gap_in_kb=max_internal_gap_in_kb,
                                                          min_roh_length=min_roh_length,
                                                          min_homozygous_snps_in_roh=min_homozygous_snps_in_roh,
                                                          min_scanning_window_hit_rate=min_scanning_window_hit_rate,
                                                          generate_overlapping_segments=generate_overlapping_segments,
                                                          max_heterozygous_snps=max_heterozygous_snps,
                                                          min_concordance_across_jointly_homozygous_variants=min_concordance_across_jointly_homozygous_variants,
                                                          homozygous_verbose=homozygous_verbose)
        self.execute(options=options)

    def test_roh_parameters(self,
                            output_dir,
                            output_prefix,
                            input_vcf_file,
                            allow_noncanonical_chromosome_names=None,
                            keep_autoconverted_files=None,
                            window_length_in_kb=None,
                            min_homozygous_snps_per_window=(1, 51, 1),
                            min_homozygous_snps_in_roh=(1, 101, 1),
                            max_heterozygous_snps_per_window=(1, 11, 1),
                            max_heterozygous_snps=(1, 21, 1),
                            max_inverse_density_of_homozygous_snps_in_kb_per_snp=(50, 1000, 50),
                            ):
        self.safe_mkdir(output_dir)
        plink_report_dict = OrderedDict()

        for i in range(*min_homozygous_snps_per_window):
            plink_report_dict[i] = OrderedDict()
            for j in range(*min_homozygous_snps_in_roh):
                plink_report_dict[j] = OrderedDict()
                for k in range(*max_heterozygous_snps_per_window):
                    plink_report_dict[k] = OrderedDict()
                    for l in range(*max_heterozygous_snps):
                        plink_report_dict[l] = OrderedDict()
                        for m in range(*max_inverse_density_of_homozygous_snps_in_kb_per_snp):

                            dir_name = "%s/%i_%i_%i_%i_%i/" % (output_dir, i, j, k, l, m)
                            description_text = "Minimum homozygous SNPs per window:\t%i\n" % i
                            description_text += "Minimum homozygous SNPs in ROh:\t%i\n" % j
                            description_text += "Maximum heterozygous SNPs per window:\t%i\n" % k
                            description_text += "Max heterozygous SNPs:\t%i\n" % l
                            description_text += "Max inverse density of homozygous SNPs(kb/SNP)" % m
                            self.safe_mkdir(dir_name, description_text=description_text, description_filename="DESCRIPTION")
                            self.find_runs_of_homozygosity("%s/%s" % (output_dir, output_prefix),
                                                           input_vcf_file=input_vcf_file,
                                                           allow_noncanonical_chromosome_names=allow_noncanonical_chromosome_names,
                                                           keep_autoconverted_files=keep_autoconverted_files,
                                                           roh_calling_method=None,
                                                           window_length_in_kb=window_length_in_kb,
                                                           min_homozygous_snps_per_window=i,
                                                           max_heterozygous_snps_per_window=k,
                                                           max_missing_snps_per_window=None,
                                                           max_inverse_density_of_homozygous_snps_in_kb_per_snp=m,
                                                           max_internal_gap_in_kb=None,
                                                           min_roh_length=None,
                                                           min_homozygous_snps_in_roh=j,
                                                           min_scanning_window_hit_rate=None,
                                                           generate_overlapping_segments=False,
                                                           max_heterozygous_snps=l,
                                                           min_concordance_across_jointly_homozygous_variants=None,
                                                           homozygous_verbose=False)
                            plink_report_dict[m] = PLINKReport("%s/%s" % (output_dir, output_prefix), report_type="ROH")

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
