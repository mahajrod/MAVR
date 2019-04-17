#!/usr/bin/env python
import os
from math import log, floor
from Tools.Abstract import Tool
from Routines import FileRoutines
from Tools.Samtools import SamtoolsV1



class STAR(Tool):
    """
    STAR doesn't understand local paths that is why os.path.abspath is used everywhere
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "STAR", path=path, max_threads=max_threads)

    # TODO: methods below are temporarly! Remove them when switch to RouToolPa will be finished
    @staticmethod
    def add_external_extraction_to_filename(filename):
        if filename[-3:] == ".gz":
            return "<(gunzip -c %s)" % filename
        elif filename[-4:] == ".bz2":
            return "<(bunzip2 -c %s)" % filename
        else:
            return filename

    @staticmethod
    def add_external_extraction_to_filelist(filelist):
        new_filelist = []

        for filename in filelist:
            if filename[-3:] == ".gz":
                new_filelist.append("<(gunzip -c %s)" % filename)
            elif filename[-4:] == ".bz2":
                new_filelist.append("<(bunzip2 -c %s)" % filename)
            else:
                new_filelist.append(filename)

        return new_filelist

    def index(self, genome_dir, genome_fasta, annotation_gtf=None, feature_from_gtf_to_use_as_exon=None, junction_tab_file=None, sjdboverhang=None,
              genomeSAindexNbases=None, genomeChrBinNbits=None, genome_size=None):

        FileRoutines.safe_mkdir(genome_dir)

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
        options += " --sjdbGTFfeatureExon %s" % feature_from_gtf_to_use_as_exon if feature_from_gtf_to_use_as_exon else ""
        self.execute(options)

    def align(self, genome_dir, forward_read_list, reverse_read_list=None, annotation_gtf=None,
              feature_from_gtf_to_use_as_exon=None, exon_tag_to_use_as_transcript_id=None,
              exon_tag_to_use_as_gene_id=None, length_of_sequences_flanking_junction=None, junction_tab_file_list=None,
              three_prime_trim=None, five_prime_trim=None, adapter_seq_for_three_prime_clip=None,
              max_mismatch_percent_for_adapter_trimming=None, three_prime_trim_after_adapter_clip=None,
              output_type="BAM", sort_bam=True, max_memory_per_thread_for_bam_sorting="4G", include_unmapped_reads_in_bam=True,
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

        forward_read_abs_path_list = [os.path.abspath(forward_read_list)] if isinstance(forward_read_list, str) else map(os.path.abspath, reverse_read_list)
        reverse_read_abs_path_list = ([os.path.abspath(reverse_read_list)] if isinstance(reverse_read_list, str) else map(os.path.abspath, reverse_read_list)) if reverse_read_list else None

        forward_read_abs_path_list = self.add_external_extraction_to_filelist(forward_read_abs_path_list)
        reverse_read_abs_path_list = self.add_external_extraction_to_filelist(reverse_read_abs_path_list)

        options += " --readFilesIn %s" % " ".join(forward_read_abs_path_list)

        options += " %s" % " ".join(reverse_read_abs_path_list) if reverse_read_abs_path_list else ""

        options += " --clip3pNbases %i" % three_prime_trim if three_prime_trim else ""
        options += " --clip5pNbases %i" % five_prime_trim if five_prime_trim else ""
        options += " --clip3pAdapterSeq %s" % adapter_seq_for_three_prime_clip if adapter_seq_for_three_prime_clip else ""
        options += " --clip3pAdapterMMp %f" % max_mismatch_percent_for_adapter_trimming if max_mismatch_percent_for_adapter_trimming else ""
        options += " --clip3pAfterAdapterNbases %i" % three_prime_trim_after_adapter_clip if three_prime_trim_after_adapter_clip else ""

        options += " --outSAMtype %s %s" % (output_type, "Unsorted") # "SortedByCoordinate" if sort_bam else "Unsorted")
        #options += " --limitBAMsortRAM %i" % max_memory_for_bam_sorting if max_memory_for_bam_sorting else ""
        options += " --outSAMunmapped Within" if include_unmapped_reads_in_bam else ""
        options += " --outReadsUnmapped Fastx" if output_unmapped_reads else ""
        options += " --outFileNamePrefix %s" % output_dir if output_dir else ""
        options += " --twopassMode Basic" if two_pass_mode else ""
        options += " --alignIntronMax %i" % max_intron_length if max_intron_length else ""

        self.execute(options)

        if sort_bam:
            print("\tSorting...")
            unsorted_bam = "%s/Aligned.out.bam" % output_dir
            sorted_bam = "%s/Aligned.sortedByCoord.out.bam" % output_dir
            SamtoolsV1.threads = self.threads
            SamtoolsV1.sort(unsorted_bam, sorted_bam, max_memory_per_thread=max_memory_per_thread_for_bam_sorting)

            print("\tIndexing bam file...")
            SamtoolsV1.index(sorted_bam)

    def align_samples(self, samples_dir, output_dir, genome_dir, genome_fasta=None, samples=None,
                      annotation_gtf=None,  sjdboverhang=None,
                      genomeSAindexNbases=None, genomeChrBinNbits=None, genome_size=None,
                      feature_from_gtf_to_use_as_exon=None, exon_tag_to_use_as_transcript_id=None,
                      exon_tag_to_use_as_gene_id=None, length_of_sequences_flanking_junction=None,
                      junction_tab_file_list=None, three_prime_trim=None, five_prime_trim=None,
                      adapter_seq_for_three_prime_clip=None, max_mismatch_percent_for_adapter_trimming=None,
                      three_prime_trim_after_adapter_clip=None, output_type="BAM", sort_bam=True,
                      max_memory_per_thread_for_bam_sorting="4G", include_unmapped_reads_in_bam=True,
                      output_unmapped_reads=True, two_pass_mode=True, max_intron_length=None):
        #STAR.threads = threads
        #STAR.path = star_dir
        
        if genome_fasta:
            STAR.index(genome_dir, genome_fasta, annotation_gtf=annotation_gtf,
                       junction_tab_file=junction_tab_file_list, sjdboverhang=sjdboverhang,
                       genomeSAindexNbases=genomeSAindexNbases, genomeChrBinNbits=genomeChrBinNbits,
                       genome_size=genome_size)
        
        sample_list = samples if samples else self.get_sample_list(samples_dir)
        
        FileRoutines.safe_mkdir(output_dir)
        
        for sample in sample_list:
            print ("Handling %s" % sample)
            sample_dir = "%s/%s/" % (samples_dir, sample)
            alignment_sample_dir = "%s/%s/" % (output_dir, sample)
            FileRoutines.safe_mkdir(alignment_sample_dir)
            filetypes, forward_files, reverse_files, se_files = FileRoutines.make_lists_forward_and_reverse_files(sample_dir)
        
            print("\tAligning reads...")
        
            self.align(genome_dir, forward_files, reverse_read_list=reverse_files,
                       annotation_gtf=annotation_gtf if not genome_fasta else None,
                       feature_from_gtf_to_use_as_exon=feature_from_gtf_to_use_as_exon,
                       exon_tag_to_use_as_transcript_id=exon_tag_to_use_as_transcript_id,
                       exon_tag_to_use_as_gene_id=exon_tag_to_use_as_gene_id,
                       length_of_sequences_flanking_junction=length_of_sequences_flanking_junction,
                       junction_tab_file_list=junction_tab_file_list,
                       three_prime_trim=three_prime_trim, five_prime_trim=five_prime_trim,
                       adapter_seq_for_three_prime_clip=adapter_seq_for_three_prime_clip,
                       max_mismatch_percent_for_adapter_trimming=max_mismatch_percent_for_adapter_trimming,
                       three_prime_trim_after_adapter_clip=three_prime_trim_after_adapter_clip,
                       output_type=output_type, sort_bam=sort_bam,
                       max_memory_per_thread_for_bam_sorting=max_memory_per_thread_for_bam_sorting,
                       include_unmapped_reads_in_bam=include_unmapped_reads_in_bam,
                       output_unmapped_reads=output_unmapped_reads, output_dir=alignment_sample_dir,
                       two_pass_mode=two_pass_mode, max_intron_length=max_intron_length)

    def align_miRNA(self, genome_dir, se_read_list, annotation_gtf=None, output_dir="./",
                    max_memory_for_bam_sorting=8000000000, max_alignments_per_read=10, no_soft_clip=True,
                    max_number_of_mismatches=1, max_relative_number_of_mismatches=None,
                    feature_name=None):
        """
        Aligns miRNA according to ENCODE pipeline
        """

        options = " --runThreadN %i" % self.threads
        options += " --genomeDir %s" % os.path.abspath(genome_dir)
        options += " --sjdbGTFfile %s" % annotation_gtf if annotation_gtf else ""

        options += " --alignEndsType EndToEnd" if no_soft_clip else ""   # End-to-end alignment, no soft-clip by default
        if max_relative_number_of_mismatches:
            options += " --outFilterMismatchNoverLmax %f" % max_relative_number_of_mismatches
        else:
            options += " --outFilterMismatchNmax %i" % max_number_of_mismatches  # Max 1 mismatch by default. maximum number of mismatches per pair, large number switches off this filter 7max number of mismatches per pair relative to read length: for 2x100b, max number of mismatches is 0.06*200=8 for the paired read
        options += " --outFilterMultimapScoreRange 0" #the score range below the maximum score for multimapping alignments
        options += " --quantMode TranscriptomeSAM GeneCounts" if annotation_gtf else ""
        options += " --outReadsUnmapped Fastx" # output of unmapped and partially mapped (i.e. mapped only one mateof a paired end read) reads in separate file(s).
        options += " --outSAMtype BAM SortedByCoordinate"
        options += " --outFilterMultimapNmax %i" % max_alignments_per_read # up to 10 alignments by default
        options += " --outSAMunmapped Within"
        options += " --outFilterScoreMinOverLread 0" # alignment will be output only if its score is higher than or equal to this valuenormalized to read length
        options += " --outFilterMatchNminOverLread 0"
        options += " --outFilterMatchNmin 16" # Output alignment if more than <value> bases were aligned
        options += " --alignSJDBoverhangMin 1000" if annotation_gtf else ""
        options += " --alignIntronMax 1" # Max intron length
        options += " --outWigType wiggle"
        options += " --outWigStrand Stranded"
        options += " --outWigNorm RPM"
        options += " --sjdbGTFfeatureExon %s" % feature_name if feature_name else ""
        options += " --outFileNamePrefix %s" % output_dir if output_dir else ""
        options += " --limitBAMsortRAM %i" % max_memory_for_bam_sorting if max_memory_for_bam_sorting else ""
        options += " --readFilesIn %s" % (os.path.abspath(se_read_list) if isinstance(se_read_list, str) else " ".join(map(os.path.abspath, se_read_list)))

        self.execute(options)
