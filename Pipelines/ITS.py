#!/usr/bin/env python

import pandas as pd

from RouToolPa.Tools.Alignment import BamUtil
from RouToolPa.Tools.Samtools import SamtoolsV1, VariantCall
from RouToolPa.Tools.Bedtools import GenomeCov
from RouToolPa.Parsers.VCF import CollectionVCF
from RouToolPa.Routines import DrawingRoutines

from Pipelines.Filtering import FilteringPipeline
from Pipelines.Alignment import AlignmentPipeline


class ITSPipeline(FilteringPipeline, AlignmentPipeline):

    def __init__(self):
        FilteringPipeline.__init__(self)
        AlignmentPipeline.__init__(self)

    def pipeline(self, samples_directory, output_directory, adapter_fragment_file, trimmomatic_adapter_file,
                 reference, index, output_prefix, aligner="bwa",
                 samples_to_handle=None, threads=4, trimmomatic_dir="", trimmer_dir="", bam_util_dir="",
                 mismatch_number=2, pe_reads_score=30, se_read_score=10,
                 min_adapter_len=1, sliding_window_size=None,
                 average_quality_threshold=15, base_quality="phred33",
                 leading_base_quality_threshold=None, trailing_base_quality_threshold=None,
                 crop_length=None, head_crop_length=None, min_len=50,
                 remove_intermediate_files=True, filtered_reads=False,
                 max_insert_size=None, max_coverage_for_variant_call=10000000, min_coverage_for_filtering=100):

        BamUtil.path = bam_util_dir
        BamUtil.threads = threads

        sample_list = samples_to_handle if samples_to_handle else self.get_sample_list(samples_directory)
        dir_dict = {"reads": {"filtered": []},
                    "alignment": {},
                    "variants": []}

        self.init_dirs(dirs=dir_dict, workdir=output_directory)
        filtered_reads_dir = "%s/reads/filtered/" % output_directory
        alignment_dir = "%s/alignment/" % output_directory

        filtered_reads_suffix = ".final"
        vcf_prefix = "%s/%s" % (output_directory, output_prefix)
        vcf_file = "%s.vcf.gz" % vcf_prefix
        tab_file = "%s.tab" % vcf_prefix

        general_stat_file = "%s/%s.filtering.stats" % (output_directory, output_prefix)

        if filtered_reads:
            filtered_reads_dir = samples_directory
        else:
            self.stirka_trimmomatic(samples_directory, filtered_reads_dir, adapter_fragment_file, trimmomatic_adapter_file,
                                    general_stat_file,
                                    samples_to_handle=sample_list, threads=threads,
                                    trimmomatic_dir=trimmomatic_dir, trimmer_dir=trimmer_dir,
                                    mismatch_number=mismatch_number, pe_reads_score=pe_reads_score,
                                    se_read_score=se_read_score,
                                    min_adapter_len=min_adapter_len, sliding_window_size=sliding_window_size,
                                    average_quality_threshold=average_quality_threshold, base_quality=base_quality,
                                    leading_base_quality_threshold=leading_base_quality_threshold,
                                    trailing_base_quality_threshold=trailing_base_quality_threshold,
                                    crop_length=crop_length, head_crop_length=head_crop_length, min_len=min_len,
                                    remove_intermediate_files=remove_intermediate_files
                                    )
        
        self.align(filtered_reads_dir, index, aligner=aligner, sample_list=sample_list,
                   outdir=alignment_dir, quality_score_type=base_quality, read_suffix=filtered_reads_suffix,
                   read_extension="fastq", alignment_format="bam",
                   threads=threads, mark_duplicates=False, platform="Illumina",
                   add_read_groups_by_picard=False,
                   gzipped_reads=False,
                   keep_inremediate_files=False,
                   mark_duplicates_tool=False,
                   calculate_coverage=True,
                   max_insert_size=max_insert_size)

        clipped_bam_list = []

        BamUtil.parallel_clipoverlap(alignment_dir, alignment_dir, sample_list, bam_suffix="", poolsize=10000000)

        for sample in sample_list:
            sample_dir = "%s/alignment/%s/" % (output_directory, sample)
            sample_prefix = "%s/%s" % (sample_dir, sample)
            raw_bam = "%s.bam" % sample_prefix
            raw_bam_coverage = "%s.tab.gz" % sample_prefix
            clipped_prefix = "%s.clipped" % sample_prefix
            clipped_bam = "%s.bam" % clipped_prefix
            clipped_bam_coverage = "%s.tab.gz" % clipped_prefix

            clipped_bam_list.append(clipped_bam)

            DrawingRoutines.draw_plot(raw_bam_coverage, sample_prefix,
                                      x_column_index=1, y_column_index=2, separator="\t",
                                      min_x=0, max_x=None, min_y=1, max_y=None, extensions=["png", ],
                                      xlabel="coverage", ylabel="position",
                                      title="Coverage of ribosomal cluster monomer by ITS lib %s" % sample,
                                      width=12, height=6,
                                      markersize=8, ylogbase=10, type="plot", grid=True,
                                      correlation=False, close_plot=True)

            GenomeCov.get_bam_coverage_stats(clipped_bam, clipped_prefix, genome_bed=None,
                                             verbose=True, calc_stats=False)
            GenomeCov.get_stats_from_coverage_file_stream_version(clipped_bam_coverage, clipped_prefix, verbose=False,
                                                                  scaffold_column=0,
                                                                  coverage_column=2,
                                                                  separator="\t",
                                                                  buffering=10000000)

            DrawingRoutines.draw_plot(clipped_bam_coverage, clipped_prefix, x_column_index=1, y_column_index=2,
                                      separator="\t",
                                      min_x=0, max_x=None, min_y=1, max_y=None, extensions=["png", ],
                                      xlabel="coverage", ylabel="position",
                                      title="Coverage of ribosomal cluster monomer by ITS lib %s" % sample, width=12,
                                      height=6,
                                      markersize=8, ylogbase=10, type="plot", grid=True,
                                      correlation=False, close_plot=True)

        VariantCall.threads = threads
        VariantCall.call_variants(reference, vcf_prefix, clipped_bam_list, chunk_length=100,
                                  split_dir="%s/split/" % output_directory,
                                  max_coverage=max_coverage_for_variant_call,
                                  min_base_quality=30, min_mapping_quality=30)

        vcf_coll = CollectionVCF(in_file=vcf_file, parsing_mode="complete")

        for sample in vcf_coll.samples:
            vcf_coll.records[(sample, "ALT_FREQ", 0)] = vcf_coll.records[sample]["AD"][1] / (
                    vcf_coll.records[sample]["AD"][1] + vcf_coll.records[sample]["AD"][0])
        idx = pd.IndexSlice
        short_coll = vcf_coll.records[vcf_coll.records[("INFO", "DP", 0)] > min_coverage_for_filtering].loc[:, idx[:, ["POS", "REF", "ALT", "ALT_FREQ"], :]]
        short_coll.columns = short_coll.columns.droplevel([1, 2])
        short_coll.index = short_coll.index.droplevel(1)
        short_coll["POS"] += 1

        short_coll.to_csv(tab_file, sep="\t", )





