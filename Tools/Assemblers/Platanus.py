#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Platanus(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "platanus", path=path, max_threads=max_threads)

    def assemble_unitigs(self, reads_files, out_prefix="out", memory_limit=None,
                         kmer_extension_step=None, branch_cutting_diference=None,
                         bubble_crush_difference=None, initial_kmer_coverage_cutoff=None,
                         min_kmer_coverage=None, kmer_extension_safety_level=None,
                         unitig_assembly_log="unitig_assembly.log"):
        # TODO: add -k and -K options

        options = " -t %i" % self.threads
        options += " -o %s" % out_prefix
        options += " -m %i" % memory_limit if memory_limit else ""
        options += " -d %f" % branch_cutting_diference if branch_cutting_diference else ""
        options += " -u %f" % bubble_crush_difference if bubble_crush_difference else ""
        options += " -s %i" % kmer_extension_step if kmer_extension_step else ""
        options += " -n %i" % initial_kmer_coverage_cutoff if initial_kmer_coverage_cutoff else ""
        options += " -c %i" % min_kmer_coverage if min_kmer_coverage else ""
        options += " -a %f" % kmer_extension_safety_level if kmer_extension_safety_level else ""
        options += " "
        options += " "
        reads_string = ""
        for reads in reads_files:
            if isinstance(reads, str):
                reads_string += " %s" % reads
            else:
                reads_string += " " + " ".join(reads)
        options += " -f %s" % reads_string
        options += " > %s" % unitig_assembly_log
        self.execute(options, cmd="platanus assemble")

    def scaffold(self, contig_file, bubble_file, pe_reads_list=None, mp_reads_list=None,
                 out_prefix="out", min_insert_size=None,scaffolding_log="scaffolding.log"):

        if (not pe_reads_list) and (not mp_reads_list):
            raise ValueError("Both PE and MP read files were not set")

        options = " -t %i" % self.threads
        options += " -o %s" % out_prefix
        options += " -c %s" % contig_file
        options += " -b %s" % bubble_file
        options += " -n %i" % min_insert_size
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " > %s" % scaffolding_log

        """
         Map paired-end (mate-pair) reads on contigs, determine the order of contigs
    and construct scaffolds.


    INPUT OPTIONS

    -ip{INT} PAIR1 [PAIR2 ...]         : Lib_id inward_pair_file (reads in 1 file, fasta or fastq)
                                         Ex. -ip1 lib1.fa

    -IP{INT} FWD1 REV1 [FWD2 REV2 ...] : Lib_id inward_pair_files (reads in 2 files, fasta or fastq)
                                         Ex. -IP1 lib1_1.fa lib1_2.fa

    -op{INT} PAIR1 [PAIR2 ...]         : Lib_id outward_pair_file (reads in 1 file, fasta or fastq)
                                         Ex. -op1 lib1.fa

    -OP{INT} FWD1 REV1 [FWD2 REV2 ...] : lib_id outward_pair_files (reads in 2 files, fasta or fastq)
                                         Ex. -OP1 lib1_1.fa lib1_2.fa

                                         The file format is automatically detected.
                                         see "***** NOTE ***** Paired-end (Mate-pair) input" below


    OTHER OPTIONS

    -a{INT1} INT2 : lib_id average_insert_size
                    Fixed average insert size (INT) is used instead of auto estimation.

    -d{INT1} INT2 : lib_id SD_insert_size
                    Fixed SD of insert size (INT) is used instead of auto estimation.

    -s INT        : Mapping seed length (default 32)
                    Seed length must not be larger than reads length. Smaller INT
                    decrease speed.

    -v INT        : Minimum overlap length (default 32)
                    If adjacent contigs have overlap (length >= INT) and properly
                    close to each other, the contigs are joined.

    -l INT        : Minimum number of link (default 3)
                    Platanus first estimates the threshold of link (number) and
                    makes scaffolds, then decreases the threshold to INT and
                    extends scaffolds.

    -u FLOAT      : Maximum difference for bubble crush (identity, default 0.1)
                    Larger FLOAT increases the number of bubbles merged. If
                    heterozygosity of the sample is high, large FLOAT may be
                    suitable (Ex. -u 0.2).


    OUTPUT FILES

    PREFIX_scaffold.fa           : assembled sequences that include gaps('N's mean gaps)

    PREFIX_scaffoldBubble.fa     : removed bubble sequences

    PREFIX_scaffoldComponent.tsv : the information about composition of scaffolds
                                   (i.e. which contigs constitute a scaffold)

--------------------------------------------------------------------------------


        """


        self.execute(options, cmd="platanus scaffold")

    def gap_close(self, tool_options):
        options = ""
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        options += " "
        self.execute(options, cmd="custom command if necessary ex: tool_command subcommand")
if __name__ == "__main__":
    pass