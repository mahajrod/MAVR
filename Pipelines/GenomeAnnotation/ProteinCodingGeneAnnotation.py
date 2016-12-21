#!/usr/bin/env python
import os
import shutil

from collections import OrderedDict

from Routines import FileRoutines
from CustomCollections.GeneralCollections import TwoLvlDict, SynDict

from Tools.Annotation import AUGUSTUS
from Tools.Filter import Cookiecutter, Trimmomatic, FaCut
from Tools.Alignment import STAR

from Parsers.FaCut import FaCutReport
from Parsers.Coockiecutter import CoockiecutterReport
from Parsers.Trimmomatic import TrimmomaticReport

from Pipelines.Filtering import FilteringPipeline


class ProteinCodingGeneAnnotation(FilteringPipeline):
    def __init__(self):
        pass

    @staticmethod
    def prepare_gene_annotation_directories(output_directory, protein_species_list=None, rnaseq_tissues_list=None):

        annotation_dir = "%s/annotation/" % output_directory

        protein_evidence_dir = "%s/protein_evidence/" % annotation_dir
        rnaseq_evidence_dir = "%s/rnaseq_evidence/" % annotation_dir
        est_evidence_dir = "%s/est_evidence/" % annotation_dir

        augustus_dir = "%s/augustus/" % annotation_dir
        augustus_hints_dir = "%s/hints/" % augustus_dir

        for directory in (annotation_dir, protein_evidence_dir, rnaseq_evidence_dir,
                          est_evidence_dir, augustus_dir, augustus_hints_dir):
            FileRoutines.save_mkdir(directory)
            #for sample in sample_list:
            #    FileRoutines.save_mkdir("%s/%s" % (directory, sample))
    """

    def star_and_htseq(self, genome_dir, samples_directory, output_directory, gff_for_htseq, count_table_file,
                       genome_fasta=None, samples_to_handle=None,
                       genome_size=None, annotation_gtf=None,
                       feature_from_gtf_to_use_as_exon=None, exon_tag_to_use_as_transcript_id=None,
                       exon_tag_to_use_as_gene_id=None, length_of_sequences_flanking_junction=None,
                       junction_tab_file_list=None,
                       three_prime_trim=None, five_prime_trim=None, adapter_seq_for_three_prime_clip=None,
                       max_mismatch_percent_for_adapter_trimming=None, three_prime_trim_after_adapter_clip=None,
                       output_type="BAM", sort_bam=True, max_memory_for_bam_sorting=None, include_unmapped_reads_in_bam=True,
                       output_unmapped_reads=True,  two_pass_mode=False, star_dir=None, threads=1, max_intron_length=None,
                       stranded_rnaseq="yes",  min_alignment_quality=10, feature_type_for_htseq="exon",
                       feature_id_attribute_for_htseq="gene_id", htseq_mode="union"):

        STAR.threads = threads
        STAR.path = star_dir

        if genome_fasta:
            STAR.index(genome_dir, genome_fasta, annotation_gtf=None, junction_tab_file=None, sjdboverhang=None,
                       genomeSAindexNbases=None, genomeChrBinNbits=None, genome_size=genome_size)

        sample_list = samples_to_handle if samples_to_handle else self.get_sample_list(samples_directory)
        self.prepare_directories(output_directory, sample_list)

        alignment_dir = "%s/alignment/" % output_directory

        count_table = TwoLvlDict()
        for sample in sample_list:
            print ("Handling %s" % sample)
            sample_dir = "%s/%s/" % (samples_directory, sample)
            alignment_sample_dir = "%s/%s/" % (alignment_dir, sample)
            filetypes, forward_files, reverse_files = FileRoutines.make_lists_forward_and_reverse_files(sample_dir)

            print "\tAligning reads..."

            STAR.align(genome_dir, forward_files, reverse_read_list=reverse_files, annotation_gtf=annotation_gtf,
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
                       max_memory_for_bam_sorting=max_memory_for_bam_sorting,
                       include_unmapped_reads_in_bam=include_unmapped_reads_in_bam,
                       output_unmapped_reads=output_unmapped_reads, output_dir=alignment_sample_dir,
                       two_pass_mode=two_pass_mode, max_intron_length=max_intron_length)

            alignment_file = "%s/Aligned.sortedByCoord.out.bam" % alignment_sample_dir

            print "\tIndexing alignment file..."
            os.system("samtools index %s" % alignment_file)

            print "\tCounting reads aligned to features..."
            count_file = "%s/%s.htseq.count" % (alignment_sample_dir, sample)

            HTSeq.count(alignment_file, gff_for_htseq, count_file, samtype="bam", order="pos",
                        stranded_rnaseq=stranded_rnaseq, min_alignment_quality=min_alignment_quality,
                        feature_type=feature_type_for_htseq, feature_id_attribute=feature_id_attribute_for_htseq,
                        mode=htseq_mode, suppress_progres_report=False)

            sample_counts = SynDict()
            sample_counts.read(count_file, header=False, separator="\t", allow_repeats_of_key=False,
                               split_values=False, values_separator=",", key_index=0, value_index=1,
                               close_after_if_file_object=False, expression=None, comments_prefix="__")
            count_table[sample] = sample_counts

        count_table.write(count_table_file)
        """
