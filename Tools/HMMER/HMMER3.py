#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil

from Bio import SearchIO

from Tools.Abstract import Tool
from Tools.LinuxTools import CGAS
from CustomCollections.GeneralCollections import IdList, SynDict
from Routines import FileRoutines


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

    @staticmethod
    def extract_hits_by_query_ids(id_list, hit_file, output_file,
                                  fileformat="tblout", close_after_if_file_object=False):
        out_fd = output_file if isinstance(output_file, file) else open(output_file, "w")
        with open(hit_file, "r") as in_fd:
            if fileformat == "domtblout" or "tblout":
                if fileformat == "domtblout":
                    query_pos_in_file = 3
                elif fileformat == "tblout":
                    query_pos_in_file = 2
                for line in in_fd:
                    if line[0] == "#":
                        out_fd.write(line)
                        continue
                    try:
                        query_id = line.split()[query_pos_in_file]
                    except IndexError:
                        print(line)
                    if query_id in id_list:
                        out_fd.write(line)

        if (not isinstance(output_file, file)) or close_after_if_file_object:
            out_fd.close()

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

    @staticmethod
    def get_ids_from_hmm3(hmmfile, ids_file=None, return_ids_list=False):
        """
        Extracts ids from hmm3 file:
            return_ids_list == True: captures output and returns ids_list, ids_file is ignored
            return_ids_list == False and ids_file == None: writes ids to stdout
            return_ids_list == False and ids_file != None: writes ids to ids_file

        """
        return CGAS.cgas(hmmfile, grep_pattern="NAME", whole_word_match=True, awk_code="{print $2}",
                         capture_output=return_ids_list, output=ids_file if not return_ids_list else None)

    def split_hmm(self, hmmfile, output_dir, num_of_recs_per_file, num_of_files=None, output_prefix=None, threads=4):

        try:
            os.mkdir(output_dir)
        except OSError:
            pass

        id_fd = CGAS.cgas(hmmfile, grep_pattern="NAME", whole_word_match=True, awk_code="{print $2}",
                          capture_output=True)

        split_index = 1
        ids_written = 0
        ids_list = IdList()
        #ids_list = read_ids(id_fd, close_after_if_file_object=False)
        ids_list.read(id_fd, close_after_if_file_object=True)
        number_of_ids = len(ids_list)
        out_prefix = self.split_filename(hmmfile)[1] if output_prefix is None else output_prefix

        num_of_ids = int(number_of_ids/num_of_files) + 1 if num_of_files else num_of_recs_per_file

        common_options = " -f"
        common_options += " %s" % hmmfile
        options_list = []
        while (ids_written + num_of_ids) <= number_of_ids:
            tmp_id_list = IdList(ids_list[ids_written:ids_written+num_of_ids])
            tmp_id_list.write("%s/%s_%i.ids" % (output_dir, out_prefix, split_index))

            options = common_options
            options += " %s/%s_%i.ids" % (output_dir, out_prefix, split_index)
            options += " > %s" % ("%s/%s_%i.hmm" % (output_dir, out_prefix, split_index))
            options_list.append(options)

            split_index += 1
            ids_written += num_of_ids

        if ids_written != number_of_ids:
            tmp_id_list = IdList(ids_list[ids_written:])
            tmp_id_list.write("%s/%s_%i.ids" % (output_dir, out_prefix, split_index))

            options = common_options
            options += " %s/%s_%i.ids" % (output_dir, out_prefix, split_index)
            options += " > %s" % ("%s/%s_%i.hmm" % (output_dir, out_prefix, split_index))
            options_list.append(options)

            split_index += 1
        #print options_list
        self.parallel_execute(options_list, cmd="hmmfetch", threads=threads)

    def __parse_hmmsxxx_common_options(self, tblout=None, domtblout=None, pfamtblout=None,
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
        #options = " --cpu %i" % self.threads
        #options += " -o %s" % outfile

        options = " --tblout %s" % tblout if tblout else ""
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

        options = self.__parse_hmmsxxx_common_options(tblout=tblout, domtblout=domtblout,
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
        options += " --cpu %i" % self.threads
        options += " -o %s" % outfile
        options += " --qformat %s" if input_format else ""
        options += " %s" % hmmfile
        options += " %s" % seqfile

        self.execute(options, cmd="hmmscan")

    def parallel_hmmscan(self, hmmfile, seqfile, output_prefix, output_dir="./",
                         num_of_seqs_per_scan=None,
                         dont_output_alignments=False,
                         model_evalue_threshold=None, model_score_threshold=None,
                         domain_evalue_threshold=None, domain_score_threshold=None,
                         model_evalue_significant_threshold=None, model_score_significant_threshold=None,
                         domain_evalue_significant_threshold=None, domain_score_significant_threshold=None,
                         use_profile_GA_gathering_cutoffs_for_thresholds=False,
                         use_profile_NC_noise_cutoffs_for_thresholds=False,
                         use_profile_TC_trusted_cutoffs_for_thresholds=False,
                         turn_off_all_heruristics=False, turn_off_bias_filter=False,
                         MSV_threshold=None, Vit_threshold=None, Fwd_threshold=None,
                         turn_off_biased_composition_score_corrections=None,
                         input_format=None, threads=None,
                         combine_output_to_single_file=True,
                         biopython_165_compartibility=False,
                         remove_tmp_dirs=True,
                         async_run=False, external_process_pool=None,
                         cpu_per_task=1,
                         handling_mode="local",
                         job_name=None,
                         log_prefix=None,
                         error_log_prefix=None,
                         max_running_jobs=None,
                         max_running_time=None,
                         cpus_per_task=None,
                         max_memmory_per_cpu=None,
                         modules_list=None,
                         environment_variables_dict=None
                         ):
        splited_fasta_dir = "%s/splited_fasta/" % output_dir
        splited_output_dir = "%s/splited_output/" % output_dir
        splited_tblout_dir = "%s/splited_tblout/" % output_dir
        splited_domtblout_dir = "%s/splited_domtblout/" % output_dir
        splited_pfamtblout_dir = "%s/splited_pfamtblout/" % output_dir

        directory_list = [
                         splited_fasta_dir,
                         splited_output_dir,
                         splited_tblout_dir,
                         splited_domtblout_dir,
                         splited_pfamtblout_dir
                         ]

        for directory in directory_list:
            self.safe_mkdir(directory)

        common_options = self.__parse_hmmsxxx_common_options(tblout=None, domtblout=None,
                                                             pfamtblout=None,
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
        common_options += " --cpu %i" % cpu_per_task
        common_options += " --qformat %s" if input_format else ""

        if handling_mode == "local":
            number_of_files = num_of_seqs_per_scan if num_of_seqs_per_scan else 5 * threads if threads else 5 * self.threads
            self.split_fasta(seqfile, splited_fasta_dir, num_of_files=number_of_files, output_prefix=output_prefix)
            input_list_of_files = sorted(os.listdir(splited_fasta_dir))
            list_of_files = []

            for filename in input_list_of_files:
                filename_prefix = self.split_filename(filename)[1]

                input_file = "%s%s" % (splited_fasta_dir, filename)
                output_file = "%s%s.hits" % (splited_output_dir, filename_prefix)
                tblout_file = "%s%s.tblout" % (splited_tblout_dir, filename_prefix) if splited_tblout_dir else None
                domtblout_file = "%s%s.domtblout" % (splited_domtblout_dir, filename_prefix) if splited_domtblout_dir else None
                pfamtblout_file = "%s%s.pfamtblout" % (splited_pfamtblout_dir, filename_prefix) if splited_pfamtblout_dir else None

                list_of_files.append((input_file, output_file, tblout_file, domtblout_file, pfamtblout_file))

            options_list = []
            out_files = []
            tblout_files = []
            domtblout_files = []
            pfamtblout_files = []

            for in_file, out_filename, tblout_file, domtblout_file, pfamtblout_file in list_of_files:
                options = common_options

                options += " --tblout %s" % tblout_file if tblout_file else ""
                options += " --domtblout %s" % domtblout_file if domtblout_file else ""
                options += " --pfamtblout %s" % pfamtblout_file if pfamtblout_file else ""
                options += " -o %s" % out_filename

                options += " %s" % hmmfile
                options += " %s" % in_file

                options_list.append(options)
                out_files.append(out_filename)
                tblout_files.append(tblout_file)
                domtblout_files.append(domtblout_file)
                pfamtblout_files.append(pfamtblout_file)

            self.parallel_execute(options_list, cmd="hmmscan", threads=threads, async_run=async_run,
                                  external_process_pool=external_process_pool)

            if combine_output_to_single_file:
                if biopython_165_compartibility:
                    CGAS.cgas(out_files, sed_string="s/^Description:.*/Description: <unknown description>/",
                              output="%s/%s.hits" % (output_dir,output_prefix))
                else:
                    CGAS.cat(out_files, output="%s/%s.hits" % (output_dir, output_prefix))
                CGAS.cat(tblout_files, output="%s/%s.tblout" % (output_dir, output_prefix))
                CGAS.cat(domtblout_files, output="%s/%s.domtblout" % (output_dir,output_prefix))
                CGAS.cat(pfamtblout_files, output="%s/%s.pfamtblout" % (output_dir,output_prefix))

            if remove_tmp_dirs:
                for tmp_dir in directory_list :
                    shutil.rmtree(tmp_dir)

        elif handling_mode == "slurm":
            number_of_files = self.split_fasta(seqfile, splited_fasta_dir, num_of_files=threads if threads else self.threads,
                                               output_prefix=output_prefix)

            slurm_cmd_options = "hmmscan %s" % common_options

            slurm_cmd_options += " --tblout %s/%s_${SLURM_ARRAY_TASK_ID}.tblout" % (splited_tblout_dir, output_prefix) if splited_tblout_dir else ""
            slurm_cmd_options += " --domtblout %s/%s_${SLURM_ARRAY_TASK_ID}.domtblout" % (splited_domtblout_dir, output_prefix) if splited_domtblout_dir else ""
            slurm_cmd_options += " --pfamtblout %s/%s_${SLURM_ARRAY_TASK_ID}.pfamtblout" % (splited_pfamtblout_dir, output_prefix) if splited_pfamtblout_dir else ""
            slurm_cmd_options += " -o %s/%s_${SLURM_ARRAY_TASK_ID}.hits" % (splited_output_dir, output_prefix)

            slurm_cmd_options += " %s" % hmmfile
            slurm_cmd_options += " %s/%s_${SLURM_ARRAY_TASK_ID}.fasta" % (splited_fasta_dir, output_prefix)

            #print number_of_files

            return self.slurm_run_job_array(job_name,
                                            log_prefix,
                                            slurm_cmd_options,
                                            error_log_prefix,
                                            "%s%s.slurm" % (output_dir, output_prefix),
                                            task_index_list=None,
                                            start_task_index=1,
                                            end_task_index=number_of_files,
                                            max_running_jobs=max_running_jobs,
                                            max_running_time=max_running_time,
                                            cpus_per_task=cpu_per_task,
                                            max_memmory_per_cpu=max_memmory_per_cpu,
                                            modules_list=modules_list,
                                            environment_variables_dict=environment_variables_dict)

            """
            self.generate_slurm_job_array_script(job_name,
                                                 log_prefix,
                                                 slurm_cmd_options,
                                                 error_log_prefix,
                                                 job_array_script_file=job_array_script_file,
                                                 task_index_list=None,
                                                 start_task_index=1,
                                                 end_task_index=number_of_files,
                                                 max_running_jobs=max_running_jobs,
                                                 max_running_time=max_running_time,
                                                 max_memmory_per_cpu=max_memmory_per_cpu)

            job_array_id = self.slurm_run_job_array(job_array_script_file)

            return job_array_id
            """

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

        options = self.__parse_hmmsxxx_common_options(tblout=tblout, domtblout=domtblout,
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
        options += " --cpu %i" % self.threads
        options += " -o %s" % outfile
        options += " -A %s" % multialignout
        options += " --tformat %s" if input_format else ""
        options += " %s" % hmmfile
        options += " %s" % seqfile

        self.execute(options, cmd="hmmsearch")

    @staticmethod
    def extract_dom_ids_hits_from_domtblout(domtblout_file, output_file=None):
        hits_dict = SynDict()
        hits_dict.read(domtblout_file, header=False, separator=None, allow_repeats_of_key=True,
                       key_index=3, value_index=1, comments_prefix="#")
        if output_file:
            hits_dict.write(output_file, splited_values=True)
        return hits_dict

    @staticmethod
    def extract_dom_names_hits_from_domtblout(domtblout_file, output_file):
        hits_dict = SynDict()
        hits_dict.read(domtblout_file, header=False, separator=None, allow_repeats_of_key=True,
                       key_index=3, value_index=0, comments_prefix="#")
        if output_file:
            hits_dict.write(output_file, splited_values=True)
        return hits_dict

    def extract_query_ids_with_hits(self, domtblout_file, output_file):
        return self.extract_ids_from_file(domtblout_file, output_file=output_file, header=False, column_separator=" ",
                                          comments_prefix="#", column_number=0)

    def extract_top_hits(self, hmmer_hits, output_prefix, parsing_mode="index_db"): #top_hits_file, top_hits_ids_file=None,not_significant_ids_file=None, not_found_ids_file=None):

        top_hits_ids = IdList()
        not_significant_ids = IdList()
        not_found_ids = IdList()

        top_hits_file = "%s.top_hits" % output_prefix
        top_hits_ids_file = "%s.top_hits.ids" % output_prefix
        not_significant_ids_file = "%s.not_significant.ids" % output_prefix
        not_found_ids_file = "%s.not_found.ids" % output_prefix

        index_file = "%s.hmmer_hits.tmp.idx" % output_prefix

        #hmm_dict = SearchIO.index_db(index_file, hmmer_hits, "hmmer3-text")

        hmm_dict = self.parse_search_file(hmmer_hits, parsing_mode, format="hmmer3-text", index_file=None)

        with open(top_hits_file, "w") as out_fd:
            out_fd.write("#query\thit\tevalue\tbitscore\n")

            for query in hmm_dict:
                if hmm_dict[query].hits:
                    if hmm_dict[query][0].is_included:
                        out_fd.write("%s\t%s\t%s\t%s\n" % (query, hmm_dict[query][0].id, hmm_dict[query][0].evalue,
                                                           hmm_dict[query][0].bitscore))
                        top_hits_ids.append(query)
                    else:
                        not_significant_ids.append(query)
                else:
                    not_found_ids.append(query)

        os.remove(index_file)

        for id_list, id_file in zip([not_significant_ids, not_found_ids_file, top_hits_ids_file],
                                    [not_significant_ids_file, not_found_ids_file, top_hits_ids_file]):
            id_file.write(id_file)


    @staticmethod
    def get_families_from_top_hits(top_hits_file, fam_file):

        hit_dict = SynDict()
        hit_dict.read(top_hits_file, allow_repeats_of_key=True, key_index=1, value_index=0, comments_prefix="#")
        hit_dict.write(fam_file, splited_values=True)

        return hit_dict


if __name__ == "__main__":
    pass
