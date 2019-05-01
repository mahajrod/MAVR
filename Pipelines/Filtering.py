#!/usr/bin/env python

import os
import shutil
from collections import OrderedDict
from RouToolPa.Tools.Filter import Cookiecutter, Trimmomatic, FaCut, Cutadapt
from RouToolPa.Parsers.FaCut import FaCutReport
from RouToolPa.Parsers.Coockiecutter import CoockiecutterReport
from RouToolPa.Parsers.Trimmomatic import TrimmomaticReport
from RouToolPa.Collections.General import TwoLvlDict
from Pipelines.Abstract import Pipeline

class FilteringPipeline(Pipeline):

    def __init__(self):
        Pipeline.__init__(self)

    def prepare_filtering_directories(self, output_directory, sample_list):
        out_dir = os.path.abspath(output_directory)
        merged_raw_dir = "%s/merged/" % out_dir
        filtered_dir = "%s/filtered/" % out_dir
        filtering_stat_dir = "%s/filtered_stat/" % out_dir
        coockie_filtered_dir = "%s/coockiecutter/" % filtered_dir
        coockie_trimmomatic_filtered_dir = "%s/coockiecutter_trimmomatic/" % filtered_dir
        coockie_trimmomatic_quality_filtered_dir = "%s/coockiecutter_trimmomatic_quality/" % filtered_dir
        final_filtered_dir = "%s/final/" % filtered_dir

        self.safe_mkdir(filtered_dir)
        for directory in merged_raw_dir, coockie_filtered_dir, coockie_trimmomatic_filtered_dir, coockie_trimmomatic_quality_filtered_dir, final_filtered_dir, filtering_stat_dir:
            self.safe_mkdir(directory)
            for sample in sample_list:
                self.safe_mkdir("%s/%s" % (directory, sample))
        return (merged_raw_dir, filtered_dir, coockie_filtered_dir, coockie_trimmomatic_filtered_dir,
                coockie_trimmomatic_quality_filtered_dir, final_filtered_dir, filtering_stat_dir)

    def filter(self, samples_directory, output_directory, adapter_fragment_file, trimmomatic_adapter_file,
               general_stat_file,
               samples_to_handle=None, threads=4, trimmomatic_dir="", coockiecutter_dir="", facut_dir="",
               mismatch_number=2, pe_reads_score=30, se_read_score=10,
               min_adapter_len=1, sliding_window_size=None,
               average_quality_threshold=15, base_quality="phred33", read_name_type="illumina",
               leading_base_quality_threshold=None, trailing_base_quality_threshold=None,
               crop_length=None, head_crop_length=None, min_len=50,
               remove_intermediate_files=False, skip_coockiecutter=False,
               retain_single_end_reads=True, input_is_se=False):

        Cookiecutter.path = coockiecutter_dir
        Trimmomatic.jar_path = trimmomatic_dir
        Trimmomatic.threads = threads
        FaCut.path = facut_dir

        self.safe_mkdir(output_directory)
        """
        merged_raw_dir = "%s/merged/" % output_directory
        filtered_dir = "%s/filtered/" % output_directory
        coockie_filtered_dir = "%s/coockiecutter/" % filtered_dir
        coockie_trimmomatic_filtered_dir = "%s/coockiecutter_trimmomatic/" % filtered_dir
        coockie_trimmomatic_quality_filtered_dir = "%s/coockiecutter_trimmomatic_quality/" % filtered_dir
        final_filtered_dir = "%s/final/" % filtered_dir
        filtering_stat_dir = "%s/filtered_stat/" % output_directory
        """
        sample_list = samples_to_handle if samples_to_handle else self.get_sample_list(samples_directory)
        merged_raw_dir, filtered_dir, coockie_filtered_dir, \
        coockie_trimmomatic_filtered_dir, coockie_trimmomatic_quality_filtered_dir, \
        final_filtered_dir, filtering_stat_dir = self.prepare_filtering_directories(output_directory, sample_list)

        filtering_statistics = TwoLvlDict()

        for sample in sample_list:
            print ("Handling sample %s" % sample)
            filtering_statistics[sample] = OrderedDict()
            merged_raw_sample_dir = "%s/%s/" % (merged_raw_dir, sample)
            #merged_forward_reads = "%s/%s_1.fq" % (merged_raw_sample_dir, sample)
            #merged_reverse_reads = "%s/%s_2.fq" % (merged_raw_sample_dir, sample)

            coockie_filtered_sample_dir = "%s/%s/" % (coockie_filtered_dir, sample)
            coockie_stats = "%s/%s.coockiecutter.stats" % (coockie_filtered_sample_dir, sample)

            coockie_trimmomatic_filtered_sample_dir = "%s/%s/" % (coockie_trimmomatic_filtered_dir, sample)

            coockie_trimmomatic_quality_filtered_sample_dir = "%s/%s/" % (coockie_trimmomatic_quality_filtered_dir, sample)
            final_filtered_sample_dir = "%s/%s/" % (final_filtered_dir, sample)
            filtering_stat_sample_dir = "%s/%s" % (filtering_stat_dir, sample)

            #"""
            print("\tMerging fastqs if necessary...")
            merged_forward_reads, merged_reverse_reads, merged_se_reads = self.combine_fastq_files(samples_directory,
                                                                                                   sample,
                                                                                                   merged_raw_sample_dir,
                                                                                                   use_links_if_merge_not_necessary=True,
                                                                                                   input_is_se=input_is_se)
            if not skip_coockiecutter:
                print("\tFiltering by Cookiecutter")
                #"""
                Cookiecutter.rm_reads(adapter_fragment_file,
                                      merged_forward_reads if merged_forward_reads else merged_se_reads,
                                      coockie_stats,
                                      right_reads=merged_reverse_reads,
                                      out_dir=coockie_filtered_sample_dir, use_dust_filter=False,
                                      dust_cutoff=None, dust_window_size=None, use_N_filter=False,
                                      read_length_cutoff=None, polyGC_length_cutoff=None)

                #"""
                print("\tParsing Cookiecutter report...")
                coockiecutter_report = CoockiecutterReport(coockie_stats, input_is_se=input_is_se)

                filtering_statistics[sample]["raw_pairs"] = coockiecutter_report.input_pairs
                filtering_statistics[sample]["pairs_after_coockiecutter"] = coockiecutter_report.retained_pairs
                filtering_statistics[sample]["pairs_after_coockiecutter,%"] = float("%.2f" % (float(coockiecutter_report.retained_pairs)/float(coockiecutter_report.input_pairs)*100))

                os.system("cp %s %s" % (coockie_stats, filtering_stat_sample_dir))

                coockie_filtered_paired_forward_reads = "%s/%s_1.ok.fastq" % (coockie_filtered_sample_dir, sample)
                coockie_filtered_paired_reverse_reads = "%s/%s_2.ok.fastq" % (coockie_filtered_sample_dir, sample)
                coockie_filtered_paired_se_reads = ""

                coockie_filtered_se_reads = "%s/%s.se.ok.fastq" % (coockie_filtered_sample_dir, sample)
            # se reads produced by Coockiecutter are ignored now!!

            #coockie_trimmomatic_filtered_sample_dir = "%s/%s/" % (coockie_trimmomatic_filtered_dir, sample)
            trimmomatic_output_prefix = "%s/%s" % (coockie_trimmomatic_filtered_sample_dir, sample)
            trimmomatic_log = "%s.trimmomatic.log" % trimmomatic_output_prefix
            #"""
            if (merged_forward_reads is None) and (merged_reverse_reads is None):
                print("Filtering by Trimmomatic...")

                Trimmomatic.filter(merged_se_reads if skip_coockiecutter else coockie_filtered_se_reads,
                                   trimmomatic_output_prefix, output_extension="fq",
                                   right_reads=None,
                                   adapters_file=trimmomatic_adapter_file,
                                   mismatch_number=mismatch_number, pe_reads_score=pe_reads_score,
                                   se_read_score=se_read_score,
                                   min_adapter_len=min_adapter_len, sliding_window_size=sliding_window_size,
                                   average_quality_threshold=average_quality_threshold,
                                   leading_base_quality_threshold=leading_base_quality_threshold,
                                   trailing_base_quality_threshold=trailing_base_quality_threshold,
                                   crop_length=crop_length, head_crop_length=head_crop_length, min_length=min_len,
                                   logfile=trimmomatic_log,
                                   base_quality=base_quality)

            else:
                print("\tFiltering by Trimmomatic...")
                Trimmomatic.filter(merged_forward_reads if skip_coockiecutter else coockie_filtered_paired_forward_reads,
                                   trimmomatic_output_prefix, output_extension="fq",
                                   right_reads=merged_reverse_reads if skip_coockiecutter else coockie_filtered_paired_reverse_reads,
                                   adapters_file=trimmomatic_adapter_file,
                                   mismatch_number=mismatch_number, pe_reads_score=pe_reads_score,
                                   se_read_score=se_read_score,
                                   min_adapter_len=min_adapter_len, sliding_window_size=sliding_window_size,
                                   average_quality_threshold=average_quality_threshold,
                                   leading_base_quality_threshold=leading_base_quality_threshold,
                                   trailing_base_quality_threshold=trailing_base_quality_threshold,
                                   crop_length=crop_length, head_crop_length=head_crop_length, min_length=min_len,
                                   logfile=trimmomatic_log,
                                   base_quality=base_quality)
            #"""
            trimmomatic_report = TrimmomaticReport(trimmomatic_log, input_is_se=input_is_se)
            if skip_coockiecutter:
                filtering_statistics[sample]["raw_pairs"] = trimmomatic_report.stats["input"]

            filtering_statistics[sample]["pairs_after_trimmomatic"] = trimmomatic_report.stats["surviving"] if input_is_se else trimmomatic_report.stats["both_surviving"]
            filtering_statistics[sample]["pairs_after_trimmomatic,%"] = trimmomatic_report.stats["surviving,%"] if input_is_se else trimmomatic_report.stats["both_surviving,%"]

            if retain_single_end_reads and not input_is_se:
                filtering_statistics[sample]["forward_se_after_trimmomatic"] = trimmomatic_report.stats["forward_only_surviving"]
                filtering_statistics[sample]["forward_se_after_trimmomatic,%"] = trimmomatic_report.stats["forward_only_surviving"]

                filtering_statistics[sample]["reverse_se_after_trimmomatic"] = trimmomatic_report.stats["reverse_only_surviving,%"]
                filtering_statistics[sample]["forward_se_after_trimmomatic,%"] = trimmomatic_report.stats["forward_only_surviving,%"]

            os.system("cp %s %s" % (trimmomatic_log, filtering_stat_sample_dir))

            coockie_trimmomatic_filtered_paired_forward_reads = "%s/%s_1.pe.fq" % (coockie_trimmomatic_filtered_sample_dir, sample)
            coockie_trimmomatic_filtered_paired_reverse_reads = "%s/%s_2.pe.fq" % (coockie_trimmomatic_filtered_sample_dir, sample)

            coockie_trimmomatic_filtered_unpaired_forward_reads = "%s/%s_1.se.fq" % (coockie_trimmomatic_filtered_sample_dir, sample)
            coockie_trimmomatic_filtered_unpaired_reverse_reads = "%s/%s_2.se.fq" % (coockie_trimmomatic_filtered_sample_dir, sample)

            coockie_trimmomatic_filtered_se_reads = "%s/%s.se.fq" % (coockie_trimmomatic_filtered_sample_dir, sample)

            final_forward_reads = "%s/%s.final_1.fastq" % (final_filtered_sample_dir, sample)
            final_reverse_reads = "%s/%s.final_2.fastq" % (final_filtered_sample_dir, sample)

            final_forward_se_reads = "%s/%s.final_1.se.fastq" % (final_filtered_sample_dir, sample)
            final_reverse_se_reads = "%s/%s.final_2.se.fastq" % (final_filtered_sample_dir, sample)

            final_se_reads = "%s/%s.final.se.fastq" % (final_filtered_sample_dir, sample)

            if sliding_window_size is None:
                facut_pe_output_prefix = "%s/%s.pe" % (coockie_trimmomatic_quality_filtered_sample_dir, sample)
                facut_forward_se_output_prefix = "%s/%s.forward.se" % (coockie_trimmomatic_quality_filtered_sample_dir, sample)
                facut_reverse_se_output_prefix = "%s/%s.reverse.se" % (coockie_trimmomatic_quality_filtered_sample_dir, sample)
                facut_pe_stat_file = "%s.facut.stat" % facut_pe_output_prefix

                facut_forward_se_stat_file = "%s.facut.stat" % facut_forward_se_output_prefix
                facut_reverse_se_stat_file = "%s.facut.stat" % facut_reverse_se_output_prefix
                #"""
                FaCut.filter_by_mean_quality(average_quality_threshold,
                                             facut_pe_output_prefix,
                                             coockie_trimmomatic_filtered_paired_forward_reads,
                                             reverse_reads=coockie_trimmomatic_filtered_paired_reverse_reads,
                                             quality_type=base_quality,
                                             stat_file=facut_pe_stat_file, name_type=read_name_type)

                FaCut.filter_by_mean_quality(average_quality_threshold,
                                             facut_forward_se_output_prefix,
                                             coockie_trimmomatic_filtered_unpaired_forward_reads,
                                             quality_type=base_quality,
                                             stat_file=facut_forward_se_stat_file, name_type=read_name_type)
                FaCut.filter_by_mean_quality(average_quality_threshold,
                                             facut_reverse_se_output_prefix,
                                             coockie_trimmomatic_filtered_unpaired_reverse_reads,
                                             quality_type=base_quality,
                                             stat_file=facut_reverse_se_stat_file, name_type=read_name_type)
                #"""
                #if input_is_se:

                #else:
                facut_report = FaCutReport(facut_pe_stat_file)

                filtering_statistics[sample]["pairs_after_facut"] = facut_report.retained_pairs
                filtering_statistics[sample]["pairs_after_facut,%"] = float("%.2f" % (float(facut_report.retained_pairs) / float(facut_report.input_pairs) * 100))
                filtering_statistics[sample]["retained_pairs_in_worst_tile,%"] = facut_report.minimum_retained_pairs_in_tiles_fraction * 100

                filtering_statistics[sample]["pairs_survived_after_filtration,%"] = float("%.2f" % (float(facut_report.retained_pairs) / filtering_statistics[sample]["raw_pairs"] * 100))

                facut_filtered_forward_reads = "%s_1.pe.fq" % facut_pe_output_prefix
                facut_filtered_reverse_reads = "%s_2.pe.fq" % facut_pe_output_prefix

                facut_filtered_forward_se_reads = "%s.se.fq" % facut_forward_se_output_prefix
                facut_filtered_reverse_se_reads = "%s.se.fq" % facut_reverse_se_output_prefix

                os.system("cp %s %s" % (facut_pe_stat_file, filtering_stat_sample_dir))
                if retain_single_end_reads:
                    os.system("cp %s %s" % (facut_forward_se_stat_file, filtering_stat_sample_dir))
                    os.system("cp %s %s" % (facut_reverse_se_stat_file, filtering_stat_sample_dir))

                os.system("ln %s %s" % (facut_filtered_forward_reads, final_forward_reads))
                os.system("ln %s %s" % (facut_filtered_reverse_reads, final_reverse_reads))
                if retain_single_end_reads and not input_is_se:
                    os.system("cat %s %s > %s" % (facut_filtered_forward_se_reads,
                                                  facut_filtered_reverse_se_reads,
                                                  final_forward_se_reads))

                    #os.system("ln %s %s" % (facut_filtered_forward_se_reads, final_forward_se_reads))
                    #os.system("ln %s %s" % (facut_filtered_reverse_se_reads, final_reverse_se_reads))

                if input_is_se:
                    pass
                    #os.system("ln %s %s" % (coockie_trimmomatic_filtered_se_reads, final_se_reads))

            else:
                os.system("ln %s %s" % (coockie_trimmomatic_filtered_paired_forward_reads, final_forward_reads))
                os.system("ln %s %s" % (coockie_trimmomatic_filtered_paired_reverse_reads, final_reverse_reads))
                if retain_single_end_reads and not input_is_se:
                    os.system("cat %s %s > %s" % (coockie_trimmomatic_filtered_unpaired_forward_reads,
                                                  coockie_trimmomatic_filtered_unpaired_reverse_reads,
                                                  final_forward_se_reads))
                    """
                    os.system("ln %s %s" % (coockie_trimmomatic_filtered_unpaired_forward_reads, final_forward_se_reads))
                    os.system("ln %s %s" % (coockie_trimmomatic_filtered_unpaired_reverse_reads, final_reverse_se_reads))

                    """
                if input_is_se:
                    os.system("ln %s %s" % (coockie_trimmomatic_filtered_se_reads, final_se_reads))
                filtering_statistics[sample]["pairs_survived_after_filtration,%"] = float("%.2f" % (float(trimmomatic_report.stats["surviving" if input_is_se else "both_surviving"]) / filtering_statistics[sample]["raw_pairs"] * 100))

            print(filtering_statistics.table_form())

            if remove_intermediate_files:
                shutil.rmtree(merged_raw_sample_dir)
                shutil.rmtree(coockie_filtered_sample_dir)
                shutil.rmtree(coockie_trimmomatic_filtered_sample_dir)
                shutil.rmtree(coockie_trimmomatic_quality_filtered_sample_dir)

        if remove_intermediate_files:
            shutil.rmtree(coockie_filtered_dir)
            shutil.rmtree(coockie_trimmomatic_filtered_dir)
            shutil.rmtree(coockie_trimmomatic_quality_filtered_dir)
            shutil.rmtree(merged_raw_dir)

        filtering_statistics.write(general_stat_file, sort=False)

    def prepare_mirna_filtering_directories(self, output_directory, sample_list):
        out_dir = os.path.abspath(output_directory)
        merged_raw_dir = "%s/merged/" % out_dir
        filtered_dir = "%s/filtered/" % out_dir
        adapter_filtered_dir = "%s/adapter_filtered/" % filtered_dir
        length_adapter_filtered_dir = "%s/length_adapter_filtered/" % filtered_dir
        final_filtered_dir = "%s/final/" % filtered_dir

        self.safe_mkdir(filtered_dir)
        for directory in merged_raw_dir, adapter_filtered_dir, length_adapter_filtered_dir, final_filtered_dir:
            self.safe_mkdir(directory)
            for sample in sample_list:
                self.safe_mkdir("%s/%s" % (directory, sample))

        return merged_raw_dir, filtered_dir, adapter_filtered_dir, length_adapter_filtered_dir, final_filtered_dir

    def filter_mirna(self, samples_directory, output_directory, adapter_list,
                     general_stat_file=None,
                     samples_to_handle=None, cutadapt_dir="",
                     min_len=15, max_len=39,
                     remove_intermediate_files=False,
                     average_quality_threshold=15, base_quality="phred33", read_name_type="illumina",
                     facut_dir="",
                     ):

        Cutadapt.path = cutadapt_dir
        FaCut.path = facut_dir

        sample_list = samples_to_handle if samples_to_handle else self.get_sample_list(samples_directory)
        merged_raw_dir, filtered_dir, adapter_filtered_dir, length_adapter_filtered_dir, final_filtered_dir = \
            self.prepare_mirna_filtering_directories(output_directory, sample_list)

        #filtering_statistics = TwoLvlDict()

        for sample in sample_list:
            print("Handling sample %s" % sample)
            #filtering_statistics[sample] = OrderedDict()
            merged_raw_sample_dir = "%s/%s/" % (merged_raw_dir, sample)

            adapter_filtered_sample_dir = "%s/%s/" % (adapter_filtered_dir, sample)
            length_adapter_filtered_sample_dir = "%s/%s/" % (length_adapter_filtered_dir, sample)
            final_filtered_sample_dir = "%s/%s/" % (final_filtered_dir, sample)


            adapter_filtered_reads_prefix = "%s/%s" % (adapter_filtered_sample_dir, sample)
            cutadapt_logfile = "%s.log" % adapter_filtered_reads_prefix
            """

            merged_forward_reads, merged_reverse_reads, merged_se_reads = self.combine_fastq_files(samples_directory,
                                                                                                   sample,
                                                                                                   merged_raw_sample_dir,
                                                                                                   use_links_if_merge_not_necessary=True)
            if merged_forward_reads or merged_reverse_reads:
                raise ValueError("Presence of PE reads in miRNAseq data!!! Stopping...")


            Cutadapt.filter(merged_se_reads, adapter_filtered_reads_prefix, reverse_file=None, format="fastq",
                            three_prime_adapter_list=adapter_list, five_prime_adapter_list=None,
                            anyway_adapter_list=None,
                            max_number_of_adapters_per_read=10, trim_Ns_on_read_end=True,
                            min_read_length_after_trimming=None,
                            logfile=cutadapt_logfile)
            """
            adapter_filtered_reads = "%s.se.fastq" % adapter_filtered_reads_prefix

            length_adapter_filtered_reads = "%s/%s.length_filtered.se.fastq" % (length_adapter_filtered_sample_dir, sample)
            length_adapter_filtered_out_reads = "%s/%s.length_filtered_out.se.fastq" % (length_adapter_filtered_sample_dir, sample)
            self.filter_se_by_length(adapter_filtered_reads,
                                     length_adapter_filtered_reads,
                                     length_adapter_filtered_out_reads,
                                     min_len=min_len, max_len=max_len)

            quality_length_adapter_filtered_prefix = "%s/%s" % (final_filtered_sample_dir, sample)
            quality_length_adapter_filtered_stat = "%s/%s.stat" % (final_filtered_sample_dir, sample)

            FaCut.filter_by_mean_quality(average_quality_threshold,
                                         quality_length_adapter_filtered_prefix,
                                         length_adapter_filtered_reads,
                                         quality_type=base_quality,
                                         stat_file=quality_length_adapter_filtered_stat, name_type=read_name_type)

            if remove_intermediate_files:
                shutil.rmtree(merged_raw_sample_dir)
                shutil.rmtree(length_adapter_filtered_sample_dir)
                shutil.rmtree(adapter_filtered_sample_dir)

        if remove_intermediate_files:
            shutil.rmtree(adapter_filtered_dir)
            shutil.rmtree(length_adapter_filtered_dir)
            shutil.rmtree(merged_raw_dir)

        #filtering_statistics.write(general_stat_file, sort=False)
