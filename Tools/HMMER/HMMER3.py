#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os

from Bio import SearchIO

from Routines.File import check_path, split_filename, read_ids
from Tools.Abstract import Tool
from Tools.LinuxTools import CGAS


class HMMER3(Tool):
    """
    http://hmmer.janelia.org/
    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "hmmbuild", path=path, max_threads=max_threads)

    def hmmbuild(self, input, output, sample_name, summary_output="hmmbuild_summary.t",
                 data_type="protein", max_insert_len=None, window_length=None,
                 alt_model_constact_strategy=None, alt_rel_seq_weigthing_strategy=None,
                 alt_eff_weighting_strategy=None, alt_prior_strategy=None, ):

        options = " --cpu %i" % self.threads
        options += " --dna" if data_type == "dna" else " --rna" if data_type == "rna" else " --amino"
        options += " -n %s" % sample_name
        options += " --maxinsertlen %i" % max_insert_len if max_insert_len else ""
        options += " --w_length %i" % window_length if window_length else ""
        options += " --%s" % alt_model_constact_strategy if alt_model_constact_strategy else ""
        options += " --%s" % alt_rel_seq_weigthing_strategy if alt_rel_seq_weigthing_strategy else ""
        options += " --%s" % alt_eff_weighting_strategy if alt_eff_weighting_strategy else ""
        options += " --%s" % alt_prior_strategy if alt_prior_strategy else ""
        options += " -o %s" % summary_output
        options += " %s" % output
        options += " %s" % input

        self.execute(options, cmd="hmmbuild")

    def hmmpress(self, hmmfile):

        options = " %s" % hmmfile

        self.execute(options, cmd="hmmpress")

    def index_hmmfile(self, hmmfile):

        options = " --index %s" % hmmfile

        self.execute(options, cmd="hmmfetch")

    def get_hmms_by_ids(self, hmmfile, id_file, output=None):

        options = " -f"
        options += " %s" % hmmfile
        options += " %s" % id_file
        options += " > %s" % output if output is not None else ""

        self.execute(options, cmd="hmmfetch")

    def get_single_hmm_by_id(self, hmmfile, hmm_id, output=None):

        options = "%s" % hmmfile
        options += " %s" % hmm_id
        options += " > %s" % output if output is not None else ""

        self.execute(options, cmd="hmmfetch")

    def hmmstat(self, hmm_file, stat_file=None):

        options = " %s" % hmm_file
        options += " > %s" % stat_file if stat_file is not None else ""

        self.execute(options, cmd="hmmstat")

    def split_hmm(self, hmmfile, output_dir, num_of_recs_per_file, num_of_files=None, output_prefix=None, threads=4):

        id_fd = CGAS.cgas(hmmfile, grep_pattern="NAME", whole_word_match=True, awk_code="{print $2}",
                                capture_output=True)

        split_index = 1
        ids_written = 0
        ids_list = read_ids(id_fd, close_after_if_file_object=True)
        number_of_ids = len(ids_list)
        out_prefix = split_filename(hmmfile)[1] if output_prefix is not None else output_prefix

        num_of_ids = int(number_of_ids/num_of_files) + 1 if num_of_files else num_of_recs_per_file

        common_options = " -f"
        common_options += " %s" % hmmfile
        options_list = []
        while (ids_written + num_of_ids) <= number_of_ids:
            options = common_options
            options += " %s" % hmmfile
            options += " > %s" % ("%s/%s_%i.hmm" % (output_dir, out_prefix, split_index))

            options_list.append(options)

            split_index += 1
            ids_written += num_of_ids

        if ids_written != number_of_ids:
            options = common_options
            options += " %s" % hmmfile
            options += " > %s" % ("%s/%s_%i.hmm" % (output_dir, out_prefix, split_index))

            options_list.append(options)

            split_index += 1

        self.parallel_execute(options_list, cmd="hmmfetch", threads=threads)

    def __parse_hmmsxxx_common_options(self, outfile, tblout=None, domtblout=None, pfamtblout=None,
                                       dont_output_alignments=False, model_evalue_threshold=None,
                                       model_score_threshold=None,
                                       domain_evalue_threshold=None, domain_score_threshold=None,
                                       model_evalue_significant_threshold=None, model_score_significant_threshold=None,
                                       domain_evalue_significant_threshold=None,
                                       domain_score_significant_threshold=None,
                                       use_profile_GA_gathering_cutoffs_for_thresholds=False,
                                       use_profile_NC_noise_cutoffs_for_thresholds=False,
                                       use_profile_TC_trusted_cutoffs_for_thresholds=False,
                                       turn_off_all_heruristics=False, turn_off_bias_filter=False,
                                       MSV_threshold=None, Vit_threshold=None, Fwd_threshold=None,
                                       turn_off_biased_composition_score_corrections=None):
        options = " --cpu %i" % self.threads
        options += " -o %s" % outfile

        options += " --tblout %s" % tblout if tblout else ""
        options += " --domtblout %s" % domtblout if domtblout else ""
        options += " --pfamtblout %s" % pfamtblout if pfamtblout else ""
        options += " --noali" if dont_output_alignments else ""
        options += " --E %f" % model_evalue_threshold if model_evalue_threshold else ""
        options += " --T %f" % model_score_threshold if model_score_threshold else ""
        options += " --domE %f" % domain_evalue_threshold if domain_evalue_threshold else ""
        options += " --domT %f" % domain_score_threshold if domain_score_threshold else ""
        options += " --incE %f" % model_evalue_significant_threshold if model_evalue_significant_threshold else ""
        options += " --incT %f" % model_score_significant_threshold if model_score_significant_threshold else ""
        options += " --domincE %f" % domain_evalue_significant_threshold if domain_evalue_significant_threshold else ""
        options += " --domincT %f" % domain_score_significant_threshold if domain_score_significant_threshold else ""
        options += " --cut_ga" if use_profile_GA_gathering_cutoffs_for_thresholds else ""
        options += " --cut_nc" if use_profile_NC_noise_cutoffs_for_thresholds else ""
        options += " --cut_tc" if use_profile_TC_trusted_cutoffs_for_thresholds else ""
        options += " --max" if turn_off_all_heruristics else ""
        options += " --nobias" if turn_off_bias_filter else ""
        options += " --F1 %f" % MSV_threshold if MSV_threshold else ""
        options += " --F2 %f" % Vit_threshold if Vit_threshold else ""
        options += " --F3 %f" % Fwd_threshold if Fwd_threshold else ""
        options += " --nonull2" if turn_off_biased_composition_score_corrections else ""
        return options

    def hmmscan(self, hmmfile, seqfile, outfile, tblout=None, domtblout=None, pfamtblout=None,
                dont_output_alignments=False, model_evalue_threshold=None, model_score_threshold=None,
                domain_evalue_threshold=None, domain_score_threshold=None,
                model_evalue_significant_threshold=None, model_score_significant_threshold=None,
                domain_evalue_significant_threshold=None, domain_score_significant_threshold=None,
                use_profile_GA_gathering_cutoffs_for_thresholds=False,
                use_profile_NC_noise_cutoffs_for_thresholds=False,
                use_profile_TC_trusted_cutoffs_for_thresholds=False,
                turn_off_all_heruristics=False, turn_off_bias_filter=False,
                MSV_threshold=None, Vit_threshold=None, Fwd_threshold=None,
                turn_off_biased_composition_score_corrections=None,
                input_format=None):

        options = self.__parse_hmmsxxx_common_options(outfile=outfile, tblout=tblout, domtblout=domtblout,
                                                      pfamtblout=pfamtblout,
                                                      dont_output_alignments=dont_output_alignments,
                                                      model_evalue_threshold=model_evalue_threshold,
                                                      model_score_threshold=model_score_threshold,
                                                      domain_evalue_threshold=domain_evalue_threshold,
                                                      domain_score_threshold=domain_score_threshold,
                                                      model_evalue_significant_threshold=model_evalue_significant_threshold,
                                                      model_score_significant_threshold=model_score_significant_threshold,
                                                      domain_evalue_significant_threshold=domain_evalue_significant_threshold,
                                                      domain_score_significant_threshold=domain_score_significant_threshold,
                                                      use_profile_GA_gathering_cutoffs_for_thresholds=use_profile_GA_gathering_cutoffs_for_thresholds,
                                                      use_profile_NC_noise_cutoffs_for_thresholds=use_profile_NC_noise_cutoffs_for_thresholds,
                                                      use_profile_TC_trusted_cutoffs_for_thresholds=use_profile_TC_trusted_cutoffs_for_thresholds,
                                                      turn_off_all_heruristics=turn_off_all_heruristics,
                                                      turn_off_bias_filter=turn_off_bias_filter,
                                                      MSV_threshold=MSV_threshold, Vit_threshold=Vit_threshold,
                                                      Fwd_threshold=Fwd_threshold,
                                                      turn_off_biased_composition_score_corrections=turn_off_biased_composition_score_corrections)

        options += " --qformat %s" if input_format else ""
        options += " %s" % hmmfile
        options += " %s" % seqfile

        self.execute(options, cmd="hmmscan")

    def parallel_hmmscan(self, hmmfile, seqfile, outfile, num_of_seqs_per_scan=None, split_dir="splited_fasta",
                         splited_output_dir="splited_output_dir",
                         tblout=None, domtblout=None, pfamtblout=None,
                         dont_output_alignments=False, model_evalue_threshold=None, model_score_threshold=None,
                         domain_evalue_threshold=None, domain_score_threshold=None,
                         model_evalue_significant_threshold=None, model_score_significant_threshold=None,
                         domain_evalue_significant_threshold=None, domain_score_significant_threshold=None,
                         use_profile_GA_gathering_cutoffs_for_thresholds=False,
                         use_profile_NC_noise_cutoffs_for_thresholds=False,
                         use_profile_TC_trusted_cutoffs_for_thresholds=False,
                         turn_off_all_heruristics=False, turn_off_bias_filter=False,
                         MSV_threshold=None, Vit_threshold=None, Fwd_threshold=None,
                         turn_off_biased_composition_score_corrections=None,
                         input_format=None):

        splited_dir = check_path(split_dir)
        number_of_files = num_of_seqs_per_scan if num_of_seqs_per_scan else 2 * self.threads
        self.split_fasta(seqfile, splited_dir, number_of_files)
        list_of_files = os.listdir(splited_dir)
        list_of_files = [("%s%s" % (splited_dir, filename), "%s%s" % (splited_output_dir, filename)) for filename in list_of_files]

        common_options = self.__parse_hmmsxxx_common_options(tblout=tblout, domtblout=domtblout,
                                                             pfamtblout=pfamtblout,
                                                             dont_output_alignments=dont_output_alignments,
                                                             model_evalue_threshold=model_evalue_threshold,
                                                             model_score_threshold=model_score_threshold,
                                                             domain_evalue_threshold=domain_evalue_threshold,
                                                             domain_score_threshold=domain_score_threshold,
                                                             model_evalue_significant_threshold=model_evalue_significant_threshold,
                                                             model_score_significant_threshold=model_score_significant_threshold,
                                                             domain_evalue_significant_threshold=domain_evalue_significant_threshold,
                                                             domain_score_significant_threshold=domain_score_significant_threshold,
                                                             use_profile_GA_gathering_cutoffs_for_thresholds=use_profile_GA_gathering_cutoffs_for_thresholds,
                                                             use_profile_NC_noise_cutoffs_for_thresholds=use_profile_NC_noise_cutoffs_for_thresholds,
                                                             use_profile_TC_trusted_cutoffs_for_thresholds=use_profile_TC_trusted_cutoffs_for_thresholds,
                                                             turn_off_all_heruristics=turn_off_all_heruristics,
                                                             turn_off_bias_filter=turn_off_bias_filter,
                                                             MSV_threshold=MSV_threshold, Vit_threshold=Vit_threshold,
                                                             Fwd_threshold=Fwd_threshold,
                                                             turn_off_biased_composition_score_corrections=turn_off_biased_composition_score_corrections)

        common_options += " --qformat %s" if input_format else ""
        options_list = []
        for in_file, out_file in list_of_files:
            options = common_options
            options += " -o %s" % out_file
            options += " %s" % hmmfile
            options += " %s" % in_file

            options_list.append(options)

        self.parallel_execute(options_list, cmd="hmmscan")

    def hmmsearch(self, hmmfile, seqfile, outfile, multialignout=None, tblout=None, domtblout=None, pfamtblout=None,
                  dont_output_alignments=False, model_evalue_threshold=None, model_score_threshold=None,
                  domain_evalue_threshold=None, domain_score_threshold=None,
                  model_evalue_significant_threshold=None, model_score_significant_threshold=None,
                  domain_evalue_significant_threshold=None, domain_score_significant_threshold=None,
                  use_profile_GA_gathering_cutoffs_for_thresholds=False,
                  use_profile_NC_noise_cutoffs_for_thresholds=False,
                  use_profile_TC_trusted_cutoffs_for_thresholds=False,
                  turn_off_all_heruristics=False, turn_off_bias_filter=False,
                  MSV_threshold=None, Vit_threshold=None, Fwd_threshold=None,
                  turn_off_biased_composition_score_corrections=None,
                  input_format=None):

        options = self.__parse_hmmsxxx_common_options(outfile=outfile, tblout=tblout, domtblout=domtblout,
                                                      pfamtblout=pfamtblout,
                                                      dont_output_alignments=dont_output_alignments,
                                                      model_evalue_threshold=model_evalue_threshold,
                                                      model_score_threshold=model_score_threshold,
                                                      domain_evalue_threshold=domain_evalue_threshold,
                                                      domain_score_threshold=domain_score_threshold,
                                                      model_evalue_significant_threshold=model_evalue_significant_threshold,
                                                      model_score_significant_threshold=model_score_significant_threshold,
                                                      domain_evalue_significant_threshold=domain_evalue_significant_threshold,
                                                      domain_score_significant_threshold=domain_score_significant_threshold,
                                                      use_profile_GA_gathering_cutoffs_for_thresholds=use_profile_GA_gathering_cutoffs_for_thresholds,
                                                      use_profile_NC_noise_cutoffs_for_thresholds=use_profile_NC_noise_cutoffs_for_thresholds,
                                                      use_profile_TC_trusted_cutoffs_for_thresholds=use_profile_TC_trusted_cutoffs_for_thresholds,
                                                      turn_off_all_heruristics=turn_off_all_heruristics,
                                                      turn_off_bias_filter=turn_off_bias_filter,
                                                      MSV_threshold=MSV_threshold, Vit_threshold=Vit_threshold,
                                                      Fwd_threshold=Fwd_threshold,
                                                      turn_off_biased_composition_score_corrections=turn_off_biased_composition_score_corrections)
        options += " -A %s" % multialignout
        options += " --tformat %s" if input_format else ""
        options += " %s" % hmmfile
        options += " %s" % seqfile

        self.execute(options, cmd="hmmsearch")


if __name__ == "__main__":
    pass