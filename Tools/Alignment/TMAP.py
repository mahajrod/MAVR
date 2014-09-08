#!/usr/bin/env python
import os
from Tools.Abstract import Tool


class TMAP(Tool):

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "tmap", path=path, max_threads=max_threads)

    def index(self, reference):

        options = ""
        options += " -f %s" % reference

        self.execute(options, cmd="tmap index")

    def parse_global_options(self, reference, reads, output=None, reads_fmt=None, output_fmt=0,
                             aln_output_mode=None, verbose=False):

        """
        -i %s               Reads formats: fa, fasta, fastq, fq, sam, bam. Default auto detection
        -o %i               Output format:
                                0: sam
                                1: compressed bam
                                2: uncompressed bam
        -a %i               Output filter for the mappings:
                                0: returns the mapping with the best score only if all other mapping havd worse
                                   score, otherwise the read is unmapped
                                1: returns mapping with the best score. If more then one mapping has this score,
                                   a random mapping with this score is returned
                                2: returns all the mappings with the best score
                                3: returns all mappings, regardless of score
        """

        options = ""
        options += " -n %i" % self.max_threads
        options += " -f %s" % reference
        options += " -r %s" % reads
        options += " -i %s" % reads_fmt if reads_fmt else ""
        options += " -s %s" % output if output else ""
        options += " -o %i" % output_fmt
        options += " -v" if verbose else ""
        options += " -a %i" % aln_output_mode if aln_output_mode else ""
        return options

    @staticmethod
    def parse_common_mapping_options(match_score=None, mismatch_penalty=None, gap_open_pen=None,
                                     gap_ext_pen=None, gap_long_pen=None, gap_long_len=None,
                                     local_band_width=None, soft_clip_type=None, duplicate_window=None,
                                     max_seed_band=None, long_hit_mult=None, score_thres=None,
                                     vsw_type=None, reads_queue_size=None,
                                     sam_read_group=None, bidirectional=False, use_seq_equal=False,
                                     unroll_banding=False, rand_read_name=False, prefix_exclude=None,
                                     suffix_exclude=None, ignore_rg_from_sam=False, input_bz2=False,
                                     input_gz=False, end_repair=None, max_adapter_bases_for_soft_clipping=None,
                                     shared_memory_key=None, min_seq_length=None, max_seq_length=None):
        # TODO: add flowspace and pairing options
        """
        -g %i               Soft clip type:
                                0:
                                1:
                                2:
                                3:

        --end-repair        End repair:
                                0: disable
                                1: prefer mismatches
                                2: prefer indels
        """
        options = ""
        options += " -A %i" % match_score if match_score else ""
        options += " -M %i" % mismatch_penalty if mismatch_penalty else ""
        options += " -O %i" % gap_open_pen if gap_open_pen else ""
        options += " -E %i" % gap_ext_pen if gap_ext_pen else ""
        options += " -G %i" % gap_long_pen if gap_long_pen else ""
        options += " -K %i" % gap_long_len if gap_long_len else ""

        options += " -w %i" % local_band_width if local_band_width else ""
        options += " -g %i" % soft_clip_type if soft_clip_type else ""
        options += " -W %i" % duplicate_window if duplicate_window else ""
        options += " -B %i" % max_seed_band if max_seed_band else ""
        options += " --long-hit-mult %f" % long_hit_mult if long_hit_mult else ""
        options += " -T %i" % score_thres if score_thres else ""
        options += " -q %i" % reads_queue_size if reads_queue_size else ""
        options += " -H %i" % vsw_type if vsw_type else ""

        options += " -U" if unroll_banding else ""
        options += " -R %s" % sam_read_group if sam_read_group else ""
        options += " -D" if bidirectional else ""
        options += " -I" if use_seq_equal else ""
        options += " -u" if rand_read_name else ""
        options += " --prefix-exclude %i" % prefix_exclude if prefix_exclude else ""
        options += " --suffix-exclude %i" % suffix_exclude if suffix_exclude else ""
        options += " -C" if ignore_rg_from_sam else ""
        options += " -j" if input_bz2 else ""
        options += " -z" if input_gz else ""
        options += " --end-repair %i" % end_repair if end_repair else ""
        options += " -J %i" % max_adapter_bases_for_soft_clipping if max_adapter_bases_for_soft_clipping else ""
        options += " -k %i" % shared_memory_key if shared_memory_key else ""
        options += " --min_seq_length %i" % min_seq_length if min_seq_length else ""
        options += " --max_seq_length %i" % max_seq_length if max_seq_length else ""

        return options

    def map1(self, reference, reads, **kwargs):
        # TODO: add parsing of  command specific options
        options = self.parse_common_mapping_options(reference, reads, **kwargs)
        self.execute(options, cmd="tmap map1")

    def map2(self, reference, reads, **kwargs):
        # TODO: add parsing of  command specific options
        options = self.parse_common_mapping_options(reference, reads, **kwargs)
        self.execute(options, cmd="tmap map2")

    def map3(self, reference, reads, **kwargs):
        # TODO: add parsing of  command specific options
        options = self.parse_common_mapping_options(reference, reads, **kwargs)
        self.execute(options, cmd="tmap map3")

    def map4(self, reference, reads, **kwargs):
        # TODO: add parsing of  command specific options
        options = self.parse_common_mapping_options(reference, reads, **kwargs)
        self.execute(options, cmd="tmap map4")

    def mapvsw(self, reference, reads, **kwargs):
        options = self.parse_common_mapping_options(reference, reads, **kwargs)
        self.execute(options, cmd="tmap mapvsw")

    def mapall(self, reference, reads, global_options_dict={}, map1_options_dict={}, map2_options_dict={}, map3_options_dict={}):
        options = "%s stage1 map1 %s  map2 %s map3 %s" % (self.parse_global_options(reference, reads,
                                                                                    **global_options_dict),
                                                          self.parse_common_mapping_options(**map1_options_dict),
                                                          self.parse_common_mapping_options(**map2_options_dict),
                                                          self.parse_common_mapping_options(**map3_options_dict))
        self.execute(options, cmd="tmap mapall")

    def fasta2pac(self):
        options = ""
        self.execute(options, cmd="tmap fasta2pac")

    def pac2bwt(self):
        options = ""
        self.execute(options, cmd="tmap pac2bwt")

    def bwt2sa(self):
        options = ""
        self.execute(options, cmd="tmap bwt2sa")

    def sff2fq(self):
        options = ""
        self.execute(options, cmd="tmap sff2fq")

    def sff2sam(self):
        options = ""
        self.execute(options, cmd="tmap sff2sam")

    def refinfo(self):
        options = ""
        self.execute(options, cmd="tmap refinfo")

    def pac2fasta(self):
        options = ""
        self.execute(options, cmd="tmap pac2fasta")

    def bwtupdate(self):
        options = ""
        self.execute(options, cmd="tmap bwtupdate")

    def indexsize(self):
        options = ""
        self.execute(options, cmd="tmap indexsize")

    def sam2fs(self):
        options = ""
        self.execute(options, cmd="tmap sam2fs")

    def sw(self):
        options = ""
        self.execute(options, cmd="tmap sw")