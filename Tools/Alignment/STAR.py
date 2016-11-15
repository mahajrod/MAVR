#!/usr/bin/env python
import os
from math import log, floor
from Tools.Abstract import Tool
from Routines import FileRoutines


class STAR(Tool):
    """
    STAR doesn't understand local paths that is why os.path.abspath is used everywhere
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "STAR", path=path, max_threads=max_threads)

    def index(self, genome_dir, genome_fasta, annotation_gtf=None, junction_tab_file=None, sjdboverhang=None,
              genomeSAindexNbases=None, genomeChrBinNbits=None, genome_size=None):

        FileRoutines.save_mkdir(genome_dir)

        options = "--runMode genomeGenerate"
        options += " --genomeDir %s" % os.path.abspath(genome_dir)
        options += " --runThreadN %i" % self.threads
        options += " --genomeFastaFiles %s" % (os.path.abspath(genome_fasta) if isinstance(genome_fasta, str) else " ".join(map(os.path.abspath, genome_fasta)))
        options += " --sjdbGTFfile %s" % annotation_gtf if annotation_gtf else ""
        options += " --sjdbFileChrStartEnd %s" % junction_tab_file if junction_tab_file else ""
        options += " --sjdbOverhang %i" % sjdboverhang if sjdboverhang else "" # number of bases taken from both sides of splice junction. 100 by default
        if genome_size:
            options += " --genomeSAindexNbases %i" % min([14, (floor(log(genome_size, 2)/2)) - 1])
        else:
            options += " --genomeSAindexNbases %i" % genomeSAindexNbases if genomeSAindexNbases else "" # size of k-mers used for preindexing of suffix array
        options += " --genomeChrBinNbits %i" % genomeChrBinNbits if genomeChrBinNbits else "" # padding size (log2) of reference sequences. 18 by default
        # recommended value min(18, log2(GenomeLength/NumberOfScaffolds))
        self.execute(options)

    def align(self, genome_dir, forward_read_list, reverse_read_list=None, annotation_gtf=None,
              feature_from_gtf_to_use_as_exon=None, exon_tag_to_use_as_transcript_id=None,
              exon_tag_to_use_as_gene_id=None, length_of_sequences_flanking_junction=None, junction_tab_file_list=None,
              three_prime_trim=None, five_prime_trim=None, adapter_seq_for_three_prime_clip=None,
              max_mismatch_percent_for_adapter_trimming=None, three_prime_trim_after_adapter_clip=None,
              output_type="BAM", sort_bam=True, max_memory_for_bam_sorting=8000000000, include_unmapped_reads_in_bam=True,
              output_unmapped_reads=True, output_dir="./", two_pass_mode=False, max_intron_length=None):
        if reverse_read_list:
            if len(forward_read_list) != len(reverse_read_list):
                raise ValueError("Wrong read file pairing")

        options = " --runThreadN %i" % self.threads
        options += " --genomeDir %s" % os.path.abspath(genome_dir)
        options += " --sjdbGTFfile %s" % annotation_gtf if annotation_gtf else ""
        options += " --sjdbGTFtagExonParentTranscript %s" % exon_tag_to_use_as_transcript_id if exon_tag_to_use_as_transcript_id else ""
        options += " --sjdbGTFtagExonParentGene %s" % exon_tag_to_use_as_gene_id if exon_tag_to_use_as_gene_id else ""
        options += " --sjdbGTFfeatureExon %s" % feature_from_gtf_to_use_as_exon if feature_from_gtf_to_use_as_exon else ""

        options += " --sjdbOverhang %i" % length_of_sequences_flanking_junction if length_of_sequences_flanking_junction else ""

        options += (" --sjdbFileChrStartEnd %s" % (os.path.abspath(junction_tab_file_list) if isinstance(junction_tab_file_list, str) else " ".join(map(os.path.abspath, junction_tab_file_list)))) if junction_tab_file_list else ""

        options += " --readFilesIn %s" % (os.path.abspath(forward_read_list) if isinstance(forward_read_list, str) else " ".join(map(os.path.abspath, forward_read_list)))
        options += " %s" % ("" if not reverse_read_list else os.path.abspath(reverse_read_list) if isinstance(reverse_read_list, str) else " ".join(map(os.path.abspath, reverse_read_list)))

        options += " --clip3pNbases %i" % three_prime_trim if three_prime_trim else ""
        options += " --clip5pNbases %i" % five_prime_trim if five_prime_trim else ""
        options += " --clip3pAdapterSeq %s" % adapter_seq_for_three_prime_clip if adapter_seq_for_three_prime_clip else ""
        options += " --clip3pAdapterMMp %f" % max_mismatch_percent_for_adapter_trimming if max_mismatch_percent_for_adapter_trimming else ""
        options += " --clip3pAfterAdapterNbases %i" % three_prime_trim_after_adapter_clip if three_prime_trim_after_adapter_clip else ""

        options += " --outSAMtype %s %s" % (output_type, "SortedByCoordinate" if sort_bam else "Unsorted")
        options += " --limitBAMsortRAM %i" % max_memory_for_bam_sorting if max_memory_for_bam_sorting else ""
        options += " --outSAMunmapped Within" if include_unmapped_reads_in_bam else ""
        options += " --outReadsUnmapped Fastx" if output_unmapped_reads else ""
        options += " --outFileNamePrefix %s" % output_dir if output_dir else ""
        options += " --twopassMode Basic" if two_pass_mode else ""
        options += " --alignIntronMax %i" % max_intron_length if max_intron_length else ""

        self.execute(options)

