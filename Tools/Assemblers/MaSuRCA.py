#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class MaSuRCA(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "masurca", path=path, max_threads=max_threads)

    def generate_config(self, libraries_list, config_file, jellyfish_hash_size=200000000, kmer_size="auto", illumina_only_assembly=True,
                        limit_jump_coverage=None, source="eukaryota", trim_long_homopolymers=False, cgwErrorRate=0.15,
                        ovlMemory="4GB", minimum_count_kmer_in_error_correction=1):
        # Library description format
        # (library_type, mean_insert_size, stdev_insert_size, read_files_list)
        # library_types: PE, SE, MP, OTHER
        config = "# DATA is specified as type {PE,JUMP,OTHER} and 5 fields:\n"
        config += "# 1)two_letter_prefix 2)mean 3)stdev 4)fastq(.gz)_fwd_reads\n"
        config += "# 5)fastq(.gz)_rev_reads. The PE reads are always assumed to be\n"
        config += "# innies, i.e. --->.<---, and JUMP are assumed to be outties\n"
        config += "# <---.--->. If there are any jump libraries that are innies, such as\n"
        config += "# longjump, specify them as JUMP and specify NEGATIVE mean. Reverse reads\n"
        config += "# are optional for PE libraries and mandatory for JUMP libraries. Any\n"
        config += "# OTHER sequence data (454, Sanger, Ion torrent, etc) must be first\n"
        config += "# converted into Celera Assembler compatible .frg files (see\n"
        config += "# http://wgs-assembler.sourceforge.com)\n"
        config += "DATA\n"

        for lib_type, mean_insert_size, stdev_insert_size, read_files_list in libraries_list:
            if lib_type == "PE":
                config += "PE= pe %i %i %s\n" % (mean_insert_size, stdev_insert_size, " ". join(read_files_list))
            elif lib_type == "MP":
                config += "JUMP= sh %i %i %s\n" % (mean_insert_size, stdev_insert_size, " ". join(read_files_list))
            else:
                for read_file in read_files_list:
                    config += "OTHER=%s\n" % read_file
        config += "END\n"
        config += "PARAMETERS\n"
        config += "#this is k-mer size for deBruijn graph values between 25 and 101 are supported, auto will compute the optimal size based on the read data and GC content\n"
        config += "GRAPH_KMER_SIZE = %s\n" % str(kmer_size)
        config += "#set this to 1 for Illumina-only assemblies and to 0 if you have 1x or more long (Sanger, 454) reads, you can also set this to 0 for large data sets with high jumping clone coverage, e.g. >50x"
        config += "USE_LINKING_MATES = %i\n" % (1 if illumina_only_assembly else 0)
        config += "#this parameter is useful if you have too many jumping library mates. Typically set it to 60 for bacteria and 300 for the other organisms"
        config += "LIMIT_JUMP_COVERAGE = %i\n" % limit_jump_coverage if limit_jump_coverage else ""
        config += "#these are the additional parameters to Celera Assembler.  do not worry about performance, number or processors or batch sizes -- these are computed automatically."
        config += "#set cgwErrorRate=0.25 for bacteria and 0.1<=cgwErrorRate<=0.15 for other organisms."
        config += "CA_PARAMETERS = cgwErrorRate=%f ovlMemory=%s\n" % (cgwErrorRate, str(ovlMemory))
        config += "#minimum count k-mers used in error correction 1 means all k-mers are used.  one can increase to 2 if coverage >100\n"
        config += "KMER_COUNT_THRESHOLD = %i\n" % minimum_count_kmer_in_error_correction
        config += "#auto-detected number of cpus to use\n"
        config += "NUM_THREADS = %i\n" % self.threads
        config += "#this is mandatory jellyfish hash size -- a safe value is estimated_genome_size*estimated_coverage\n"
        config += "JF_SIZE = %i\n" % jellyfish_hash_size
        config += "#this specifies if we do (1) or do not (0) want to trim long runs of homopolymers (e.g. GGGGGGGG) from 3' read ends, use it for high GC genomes\n"
        config += "DO_HOMOPOLYMER_TRIM = %i\n" % (1 if trim_long_homopolymers else 0)
        config += "END\n"

        with open(config_file, "w") as out_fd:
            out_fd.write(config)

    def generate_script(self, config_file, output_script):

        options = " %s" % config_file
        options += " -o %s" % output_script

        self.execute(options)

    def assembly(self, libraries_list, config_file, output_script, jellyfish_hash_size=200000000, kmer_size="auto", illumina_only_assembly=True,
                        limit_jump_coverage=None, source="eukaryota", trim_long_homopolymers=False, cgwErrorRate=0.15,
                        ovlMemory="4GB", minimum_count_kmer_in_error_correction=1):

        self.generate_config(libraries_list=libraries_list, config_file=config_file,
                             jellyfish_hash_size=jellyfish_hash_size, kmer_size=kmer_size,
                             illumina_only_assembly=illumina_only_assembly, limit_jump_coverage=limit_jump_coverage,
                             source=source, trim_long_homopolymers=trim_long_homopolymers, cgwErrorRate=cgwErrorRate,
                             ovlMemory=ovlMemory, minimum_count_kmer_in_error_correction=minimum_count_kmer_in_error_correction)
        self.generate_script(config_file, output_script)
        self.execute("", cmd=output_script)


if __name__ == "__main__":
    pass