#!/usr/bin/env python
import os
from RouToolPa.Tools.Alignment import STAR
from RouToolPa.Tools.Expression import HTSeq
from RouToolPa.Collections.General import TwoLvlDict, SynDict
from Pipelines.Filtering import FilteringPipeline



class DiffExpressionPipeline(FilteringPipeline):

    def __init__(self):
        FilteringPipeline.__init__(self)

    def prepare_diff_expression_directories(self, output_directory, sample_list):

        alignment_dir = "%s/alignment/" % output_directory

        for directory in (alignment_dir, ):
            self.safe_mkdir(directory)
            for sample in sample_list:
                self.safe_mkdir("%s/%s" % (directory, sample))

    def star_and_htseq(self, genome_dir, samples_directory, output_directory, gff_for_htseq, count_table_file_prefix,
                       genome_fasta=None, samples_to_handle=None,
                       genome_size=None, annotation_gtf=None,
                       feature_from_gtf_to_use_as_exon=None, exon_tag_to_use_as_transcript_id=None,
                       exon_tag_to_use_as_gene_id=None, length_of_sequences_flanking_junction=None,
                       junction_tab_file_list=None,
                       three_prime_trim=None, five_prime_trim=None, adapter_seq_for_three_prime_clip=None,
                       max_mismatch_percent_for_adapter_trimming=None, three_prime_trim_after_adapter_clip=None,
                       output_type="BAM", sort_bam=True, max_memory_per_thread_for_bam_sorting="4G",
                       include_unmapped_reads_in_bam=True, output_unmapped_reads=True,
                       two_pass_mode=False, star_dir=None, threads=1, max_intron_length=None,
                       stranded_rnaseq="yes",  min_alignment_quality=10, feature_type_for_htseq="exon",
                       feature_id_attribute_for_htseq="gene_id", htseq_mode="union"):

        STAR.threads = threads
        STAR.path = star_dir

        if genome_fasta:
            STAR.index(genome_dir, genome_fasta, annotation_gtf=None, junction_tab_file=None, sjdboverhang=None,
                       genomeSAindexNbases=None, genomeChrBinNbits=None, genome_size=genome_size)

        sample_list = samples_to_handle if samples_to_handle else self.get_sample_list(samples_directory)
        self.prepare_diff_expression_directories(output_directory, sample_list)

        alignment_dir = "%s/alignment/" % output_directory

        count_pe_table = TwoLvlDict()
        count_se_table = TwoLvlDict()
        count_all_table = TwoLvlDict()
        count_pe_table_file = "%s/%s.pe.tab" % (output_directory, count_table_file_prefix)
        count_se_table_file = "%%s/%s.se.tab" % (output_directory, count_table_file_prefix)
        count_all_table_file = "%s/%s.all.tab" % (output_directory, count_table_file_prefix)

        for sample in sample_list:
            print ("Handling %s" % sample)
            sample_dir = "%s/%s/" % (samples_directory, sample)
            alignment_sample_dir = "%s/%s/" % (alignment_dir, sample)
            alignment_sample_se_dir = "%s/se/" % alignment_sample_dir
            filetypes, forward_files, reverse_files, se_files = self.make_lists_forward_and_reverse_files(sample_dir)

            if se_files:
                self.safe_mkdir(alignment_sample_se_dir)

            print("\tAligning paired reads...")
            count_file = "%s/%s.htseq.count" % (alignment_sample_dir, sample)
            #"""
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
                       max_memory_per_thread_for_bam_sorting=max_memory_per_thread_for_bam_sorting,
                       include_unmapped_reads_in_bam=include_unmapped_reads_in_bam,
                       output_unmapped_reads=output_unmapped_reads, output_dir=alignment_sample_dir,
                       two_pass_mode=two_pass_mode, max_intron_length=max_intron_length)

            alignment_file = "%s/Aligned.sortedByCoord.out.bam" % alignment_sample_dir

            print("\tIndexing alignment file for paired reads...")
            os.system("samtools index %s" % alignment_file)

            print("\tCounting paired reads aligned to features...")


            HTSeq.count(alignment_file, gff_for_htseq, count_file, samtype="bam", order="pos",
                        stranded_rnaseq=stranded_rnaseq, min_alignment_quality=min_alignment_quality,
                        feature_type=feature_type_for_htseq, feature_id_attribute=feature_id_attribute_for_htseq,
                        mode=htseq_mode, suppress_progres_report=False)
            #"""
            sample_counts = SynDict(filename=count_file, header=False, separator="\t", allow_repeats_of_key=False,
                                    split_values=False, values_separator=",", key_index=0, value_index=1,
                                    close_after_if_file_object=False, expression=int, comments_prefix="__")
            count_pe_table[sample] = sample_counts

            if se_files:
                print("\tAligning single reads...")
                count_se_file = "%s/%s.htseq.count" % (alignment_sample_se_dir, sample)
                #"""
                STAR.align(genome_dir, se_files, reverse_read_list=None, annotation_gtf=annotation_gtf,
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
                           output_unmapped_reads=output_unmapped_reads, output_dir=alignment_sample_se_dir,
                           two_pass_mode=two_pass_mode, max_intron_length=max_intron_length)

                alignment_se_file = "%s/Aligned.sortedByCoord.out.bam" % alignment_sample_se_dir

                print("\tIndexing alignment file for single reads...")
                os.system("samtools index %s" % alignment_se_file)

                print("\tCounting single reads aligned to features...")

                HTSeq.count(alignment_se_file, gff_for_htseq, count_se_file, samtype="bam", order="pos",
                            stranded_rnaseq=stranded_rnaseq, min_alignment_quality=min_alignment_quality,
                            feature_type=feature_type_for_htseq, feature_id_attribute=feature_id_attribute_for_htseq,
                            mode=htseq_mode, suppress_progres_report=False)
                #"""

                sample_se_counts = SynDict(filename=count_se_file, header=False, separator="\t", allow_repeats_of_key=False,
                                           split_values=False, values_separator=",", key_index=0, value_index=1,
                                           close_after_if_file_object=False, expression=int, comments_prefix="__")

                count_se_table[sample] = sample_se_counts
            else:
                count_se_table[sample] = SynDict()
            count_all_table[sample] = SynDict()
            if se_files:
                for gene_id in set(sample_counts.keys()) | set(sample_se_counts.keys()):
                    if (gene_id in sample_counts) and (gene_id in sample_se_counts):
                        count_all_table[sample][gene_id] = sample_counts[gene_id] + sample_se_counts[gene_id]
                    elif gene_id in sample_counts:
                        count_all_table[sample][gene_id] = sample_counts[gene_id]
                    elif gene_id in sample_se_counts:
                        count_all_table[sample][gene_id] = sample_se_counts[gene_id]
            else:
                count_all_table[sample] = count_pe_table[sample]

        count_pe_table.write(count_pe_table_file)
        count_se_table.write(count_se_table_file)
        count_all_table.write(count_all_table_file)
