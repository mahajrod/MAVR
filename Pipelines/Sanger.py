#!/usr/bin/env python

import os
from collections import OrderedDict

import numpy as np

from Bio import SeqIO

from Pipelines.Abstract import Pipeline
from RouToolPa.Collections.General import IdList, TwoLvlDict
from RouToolPa.Routines import MatplotlibRoutines

class SangerPipeline(Pipeline):

    def __init__(self, workdir="./", max_threads=1):
        Pipeline.__init__(self, workdir=workdir, max_threads=max_threads)

        self.dirs = {"fastq": {
                               "raw": [],
                               "trimmed": [],
                               },
                     "fasta": {
                               "raw": [],
                               "trimmed": [],
                              },
                     "qual_plot": {
                                   "raw": [],
                                   "trimmed": [],
                                  },
                     }

        self.sanger_extention_list = [".ab1"]

    @staticmethod
    def is_sanger_file(filename):
        if not os.path.isdir(filename) and ((filename[-4:] == ".ab1") or (filename[-7:] == ".ab1.gz") or (filename[-7:] == ".ab1.bz2")):
            return True
        return False

    def handle_sanger_data(self, input_dir, output_prefix, outdir=None, read_subfolders=False,
                           min_mean_qual=0, min_median_qual=0, min_len=50):
        if outdir:
            self.workdir = outdir

        self.init_dirs()

        sanger_filelist = self.make_list_of_path_to_files(input_dir,
                                                          expression=self.is_sanger_file,
                                                          recursive=read_subfolders,
                                                          return_absolute_paths=True)
        stat_dict = TwoLvlDict()
        record_dict = OrderedDict()
        trimmed_record_dict = OrderedDict()
        excluded_list = IdList()
        excluded_counter = 0
        low_quality_counter = 0
        too_short_counter = 0

        merged_raw_fastq = "%s/%s.raw.fastq" % (self.workdir, output_prefix)
        merged_raw_fasta = "%s/%s.raw.fasta" % (self.workdir, output_prefix)
        merged_trimmed_fastq = "%s/%s.trimmed.fastq" % (self.workdir, output_prefix)
        merged_trimmed_fasta = "%s/%s.trimmed.fasta" % (self.workdir, output_prefix)

        for filename in sanger_filelist:
            filename_list = self.split_filename(filename)

            record_raw_fastq = "%s/fastq/raw/%s.raw.fastq" % (self.workdir, filename_list[1])
            record_raw_fasta = "%s/fasta/raw/%s.raw.fasta" % (self.workdir, filename_list[1])
            record_raw_qual_plot_prefix = "%s/qual_plot/raw/%s.raw.qual" % (self.workdir, filename_list[1])

            record_trimmed_fastq = "%s/fastq/trimmed/%s.trimmed.fastq" % (self.workdir, filename_list[1])
            record_trimmed_fasta = "%s/fasta/trimmed/%s.trimmed.fasta" % (self.workdir, filename_list[1])
            record_trimmed_qual_plot_prefix = "%s/qual_plot/trimmed/%s.trimmed.qual" % (self.workdir, filename_list[1])

            record = SeqIO.read(self.metaopen(filename, "rb"), format="abi")
            record_dict[record.id] = record
            SeqIO.write(record, record_raw_fastq, format="fastq")
            SeqIO.write(record, record_raw_fasta, format="fasta")

            trimmed_record = SeqIO.AbiIO._abi_trim(record)

            stat_dict[record.id] = OrderedDict({
                                                "raw_len": len(record),
                                                "raw_mean_qual": np.mean(record.letter_annotations["phred_quality"]),
                                                "raw_median_qual": np.median(record.letter_annotations["phred_quality"]),
                                                "trimmed_len": len(trimmed_record),
                                                "trimmed_mean_qual": np.mean(trimmed_record.letter_annotations["phred_quality"]),
                                                "trimmed_median_qual": np.median(trimmed_record.letter_annotations["phred_quality"]),
                                                "retained": "-",
                                                })
            MatplotlibRoutines.draw_bar_plot(record.letter_annotations["phred_quality"], record_raw_qual_plot_prefix,
                                             extentions=["png"], xlabel="Position", ylabel="Phred quality",
                                             title="Per base quality", min_value=None, max_value=None, new_figure=True,
                                             figsize=(3 * (int(len(record) / 100) + 1), 3), close_figure=True)

            if stat_dict[record.id]["trimmed_len"] >= min_len:
                if min_median_qual:
                    if (stat_dict[record.id]["trimmed_median_qual"] >= min_median_qual) and (stat_dict[record.id]["trimmed_mean_qual"] >= min_mean_qual):
                        stat_dict[record.id]["retained"] = "+"
                    else:
                        low_quality_counter += 1
                else:
                    stat_dict[record.id]["retained"] = "+"
            else:
                too_short_counter += 1

            if stat_dict[record.id]["retained"] == "-":
                excluded_list.append(record.id)
                continue

            SeqIO.write(trimmed_record, record_trimmed_fastq, format="fastq")
            SeqIO.write(trimmed_record, record_trimmed_fasta, format="fasta")

            MatplotlibRoutines.draw_bar_plot(trimmed_record.letter_annotations["phred_quality"],
                                             record_trimmed_qual_plot_prefix,
                                             extentions=["png"], xlabel="Position", ylabel="Phred quality",
                                             title="Per base quality", min_value=None, max_value=None, new_figure=True,
                                             figsize=(3 * (int(len(record) / 100) + 1), 3),
                                             close_figure=True)

            trimmed_record_dict[record.id] = trimmed_record

        SeqIO.write(self.record_from_dict_generator(record_dict), merged_raw_fastq, format="fastq")
        SeqIO.write(self.record_from_dict_generator(record_dict), merged_raw_fasta, format="fasta")

        SeqIO.write(self.record_from_dict_generator(trimmed_record_dict), merged_trimmed_fastq, format="fastq")
        SeqIO.write(self.record_from_dict_generator(trimmed_record_dict), merged_trimmed_fasta, format="fasta")

        excluded_list.write("%s.excluded.ids" % output_prefix)
        stat_dict.write(out_filename="%s.stats" % output_prefix)

        print("Excluded: %i" % excluded_counter)
        print("\tToo short( < %i ): %i" % (min_len, too_short_counter))
        print("\tLow quality( median < %i or mean < %i ): %i" % (min_median_qual, min_mean_qual, low_quality_counter))
