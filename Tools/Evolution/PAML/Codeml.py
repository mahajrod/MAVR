#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import sys

import multiprocessing as mp

from collections import OrderedDict

from scipy.stats import chisqprob
from statsmodels.sandbox.stats.multicomp import fdrcorrection0

from Tools.Abstract import Tool
from Parsers.PAML import CodeMLReport
from Routines import FileRoutines

from CustomCollections.GeneralCollections import TwoLvlDict

print_mutex = mp.Lock()


def extract_trees_from_codeml_report(list_of_options):
    # this function is global because of stutid damned pickle mode in python!!!!!
    # use mutex to safe write to stdout from multiple threads
    # list of options = [sample_directory, tree_file, report_suffix, output_prefix]
    #sys.stdout.write("Handling %s\n" % list_of_options[0])
    #print list_of_options
    sample_name = list_of_options[0].split("/")
    sample_name = sample_name[-1] if sample_name[-1] != "" else sample_name[-2]
    work_dir = os.getcwd()
    sample_dir = os.path.abspath(list_of_options[0])
    #print sample_dir
    list_of_files = os.listdir(sample_dir)
    report_files_list = []
    suffix_length = len(list_of_options[2])
    for filename in list_of_files:
        if len(filename) < suffix_length:
            continue
        if filename[-suffix_length:] == list_of_options[2]:
            report_files_list.append(filename)

    print_string = "Handling %s\n" % list_of_options[0]
    #print "XXXX"

    os.chdir(sample_dir)
    codeml_report = CodeMLReport(report_files_list[0], treefile=list_of_options[1])

    error_code = 0
    if not report_files_list:
        print_string += "\tNo report file with suffix %s were found\n\tSkipping\n" % list_of_options[2]
        error_code = -1
    elif len(report_files_list) > 1:
        print_string += "\tWere found several (%i) report files with suffix %s\n\tSkipping\n" % (len(report_files_list),
                                                                                                 list_of_options[2])
        error_code = -2
    elif (codeml_report.dNtree is None) or (codeml_report.dStree is None) or (codeml_report.Wtree is None):
        print_string += "\tCodeML report is not full\n\tSkipping\n"
        error_code = -3

    print_mutex.acquire()
    sys.stdout.write(print_string)
    print_mutex.release()
    if error_code < 0:
        return sample_name, error_code

    codeml_report.write_trees(list_of_options[3])
    codeml_report.get_all_values(list_of_options[3] + ".all.values")
    codeml_report.get_feature_values(mode="leaves")
    os.chdir(work_dir)

    tmp = codeml_report.find_leaves_with_positive_selection()
    #if tmp:
    #    print sample_name, tmp
    extract_trees_from_codeml_report.queue.put((sample_name, codeml_report.find_leaves_with_positive_selection()))
    return sample_name, codeml_report.find_leaves_with_positive_selection()


def extract_trees_from_codeml_report_init(queue):
    extract_trees_from_codeml_report.queue = queue


def results_extraction_listener(queue, output_file_prefix, selected_species_list=None):
    """listens for messages on the queue, writes to file."""

    positive_selection_dict = TwoLvlDict()
    selected_species_positive_selection_dict = TwoLvlDict()
    error_fd = open("errors.err", "w")
    error_fd.write("#sample\terror_code\n")
    while 1:
        result = queue.get()
        if isinstance(result[1], int):
            error_fd.write("%s\t%i\n" % (result[0], result[1]))
            continue
        if result == 'finish':
            print "AAA"
            positive_selection_dict.write("%s.all" % output_file_prefix, absent_symbol=".")
            if selected_species_list:
                selected_species_positive_selection_dict.write("%s.selected_species" % output_file_prefix,
                                                               absent_symbol=".")
            # print positive_selection_dict.table_form(absent_symbol=".")
            break
        if result[1]:
            positive_selection_dict[result[0]] = result[1]
            if selected_species_list:
                for species in selected_species_list:
                    if species in result[1]:
                        print species
                        if result[0] not in selected_species_positive_selection_dict:
                            selected_species_positive_selection_dict[result[0]] = {}
                        selected_species_positive_selection_dict[result[0]][species] = result[1][species]


class Codeml(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "codeml", path=path, max_threads=max_threads)
        self.example = """
        treefile = ../tree_no_support.nwk
        seqfile = ../merged_alignment.phy
        outfile = model_0           * main result file name
        noisy = 3  * 0,1,2,3,9: how much rubbish on the screen
        verbose = 1  * 0: concise; 1: detailed, 2: too much
        runmode = 0

        seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
        CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table
        clock = 0
        aaDist = 1  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
        model = 0
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)
        NSsites = 0
        icode = 1  * 0:universal code; 1:mammalian mt; 2-10:see below
        fix_kappa = 0
        kappa = 5
        fix_omega = 0
        omega = 0.2

        getSE = 0
        RateAncestor = 0

        Small_Diff = .5e-5
        cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
        method = 0   * 0: simultaneous; 1: one branch at a time
        """

    def write_example_ctl(self, example_file):
        with open(example_file, "w") as example_fd:
            example_fd.write(self.example)

    def print_example_ctl(self):
        print(self.example)

    @staticmethod
    def generate_ctl_file(seq_file, tree_file, out_file, ctl_file, seq_type="codons", codon_frequency="F3X4", noisy=3,
                          verbose="concise", runmode=0, clock=0, aminoacid_distance=None, model=1, nssites=0,
                          genetic_code=0, fix_kappa=False, kappa=5, fix_omega=False, omega=0.2, fix_alpha=True, alpha=0,
                          getSE=0, RateAncestor=0, small_difference=0.00001, clean_data=True, method=0, Mgene=None):

        if verbose == "concise":
            verb = 0
        elif verbose == "detailed":
            verb = 1
        elif verbose == "much":
            verb = 2
        else:
            verb = 0

        if seq_type == "codons":
            seq_t = 1
        elif seq_type == "aminoacids":
            seq_t = 2
        elif seq_type == "translate":
            seq_t = 3
        else:
            raise ValueError("Unknown type of sequence was set.")

        if codon_frequency == "equal":
            codon_freq = 0
        elif codon_frequency == "F1X4":
            codon_freq = 1
        elif codon_frequency == "F3X4":
            codon_freq = 2
        elif codon_frequency == "table":
            codon_freq = 3
        else:
            raise ValueError("Codon frequency was not recognized")

        options = "treefile = %s\n" % tree_file
        options += "seqfile = %s\n" % seq_file
        options += "outfile = %s\n" % out_file
        options += "noisy = %i  * 0,1,2,3,9: how much rubbish on the screen\n" % noisy
        options += "verbose = %i  * 0: concise; 1: detailed, 2: too much\n" % verb
        options += "runmode = %i * 0: for this mode only unrooted tree is allowed; 1,2: only rooted tree\n\n" % runmode

        options += "seqtype = %i  * 1:codons; 2:AAs; 3:codons-->AAs\n" % seq_t
        options += "CodonFreq = %i  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table\n" % codon_freq
        options += "clock = %s\n" % clock
        options += "aaDist = %i  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a\n" % aminoacid_distance if aminoacid_distance else ""
        options += """model = %i
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)\n""" % model
        options += "NSsites = %i\n" % nssites
        options += "icode = %i  * 0:universal code; 1:mammalian mt; 2-10:see below\n" % genetic_code
        options += "fix_kappa = %i\n" % (1 if fix_kappa else 0)
        options += "kappa = %i\n" % kappa
        options += "fix_omega = %i\n" % (1 if fix_omega else 0)
        options += "omega = %f\n\n" % omega

        options += "fix_alpha = %i\n" % (1 if fix_alpha else 0)
        options += "alpha = %f\n" % alpha


        options += "getSE = %i\n" % getSE
        options += "RateAncestor = %i\n" % RateAncestor

        options += "Small_Diff = %f\n" % small_difference
        options += "cleandata = %i  * remove sites with ambiguity data (1:yes, 0:no)?\n" % (1 if clean_data else 0)
        options += "method = %i   * 0: simultaneous; 1: one branch at a time\n" % method
        options += ("""Mgene = %i * is used only with G option in sequence file.
                   * 0: same k and pi, different cs(proportional branch lengths between genes)
                   * 1: different k, pi and unproportional branch length
                   * 2: same k, different pi and cs
                   * 3: same pi, different k and cs
                   * 4: different k, pi, cs(proportional branch length)\n""" % Mgene) if Mgene else ""

        with open(ctl_file, "w") as ctl_fd:
            ctl_fd.write(options)

    def run_codeml(self, ctl_file):

        options = " %s" % ctl_file

        self.execute(options=options)

    def parallel_codeml(self, in_dir, tree_file, out_dir, seq_type="codons", codon_frequency="F3X4", noisy=0,
                          verbose="concise", runmode=0, clock=0, aminoacid_distance=None, model=1, nssites=0,
                          genetic_code=0, fix_kappa=False, kappa=5, fix_omega=False, omega=0.2, getSE=0, RateAncestor=0,
                          small_difference=0.000001, clean_data=True, method=0, Mgene=None):

        FileRoutines.safe_mkdir(out_dir)
        alignment_files_list = FileRoutines.make_list_of_path_to_files(in_dir)
        tree_file_abs_path = os.path.abspath(tree_file)
        options_list = []
        dir_list = []
        for filename in alignment_files_list:
            directory, basename, extension = FileRoutines.split_filename(filename)
            filename_out_dir = os.path.abspath("%s/%s/" % (out_dir, basename))
            out_file = "%s/%s.out" % (filename_out_dir, basename)
            ctl_file = "%s/%s.ctl" % (filename_out_dir, basename)

            options_list.append(ctl_file)
            dir_list.append(filename_out_dir)
            FileRoutines.safe_mkdir(filename_out_dir)
            self.generate_ctl_file(os.path.abspath(filename), tree_file_abs_path, out_file, ctl_file, seq_type=seq_type,
                                   codon_frequency=codon_frequency, noisy=noisy, verbose=verbose, runmode=runmode,
                                   clock=clock, aminoacid_distance=aminoacid_distance, model=model, nssites=nssites,
                                   genetic_code=genetic_code, fix_kappa=fix_kappa, kappa=kappa, fix_omega=fix_omega,
                                   omega=omega, getSE=getSE, RateAncestor=RateAncestor, Mgene=Mgene,
                                   small_difference=small_difference, clean_data=clean_data, method=method)
        self.parallel_execute(options_list, dir_list=dir_list)

    def parallel_results_extraction(self, in_dir, tree_file, report_suffix, output_prefix,
                                    report_file_prefix, selected_species_list=None):
        work_dir = os.getcwd()
        input_directory = os.path.abspath(in_dir)
        dir_list = sorted(os.listdir(input_directory))
        tree_file_abs_path = os.path.abspath(tree_file)
        options_list = []

        for filename in dir_list:
            sample_dir = "%s/%s/" % (input_directory, filename)
            if os.path.isdir(sample_dir):
                options_list.append([sample_dir, tree_file_abs_path, report_suffix, output_prefix])

        os.chdir(input_directory)

        manager = mp.Manager()
        queue = manager.Queue()
        process_pool = mp.Pool(self.threads + 1, extract_trees_from_codeml_report_init, [queue])
        watcher = process_pool.apply_async(results_extraction_listener,
                                           (queue, report_file_prefix, selected_species_list))

        #print options_list
        process_pool.map(extract_trees_from_codeml_report, options_list)
        queue.put('finish')
        process_pool.close()
        os.chdir(work_dir)

    def parallel_positive_selection_test(self, in_dir, tree_file, out_dir, results_file, seq_type="codons", codon_frequency="F3X4", noisy=3,
                                         verbose="concise", runmode=0, clock=0, aminoacid_distance=None,
                                         genetic_code=0, fix_kappa=False, kappa=5,  getSE=0, RateAncestor=0,
                                         small_difference=0.000001, clean_data=True, method=0):
        """
        This function implements positive selection test (branch-site model)
        for branch labeled in tree file using model_A vs model_A_null(omega fixed to 1) comparison
        """

        FileRoutines.safe_mkdir(out_dir)
        alignment_files_list = FileRoutines.make_list_of_path_to_files(in_dir)
        tree_file_abs_path = os.path.abspath(tree_file)
        options_list = []
        dir_list = []
        basename_dir_list = []
        model_list = ["Model_A", "Model_A_null"]
        fix_omega_dict = {"Model_A": False, "Model_A_null": True}
        for filename in alignment_files_list:
            directory, basename, extension = FileRoutines.split_filename(filename)
            filename_out_dir = os.path.abspath("%s/%s/" % (out_dir, basename))
            basename_dir_list.append(basename)
            FileRoutines.safe_mkdir(filename_out_dir)

            for model in model_list:
                model_dir = "%s/%s/" % (filename_out_dir, model)
                FileRoutines.safe_mkdir(model_dir)
                out_file = "%s/%s/%s.out" % (filename_out_dir, model, basename)
                ctl_file = "%s/%s/%s.ctl" % (filename_out_dir, model, basename)

                options_list.append("%s.ctl" % basename)
                dir_list.append(model_dir)

                self.generate_ctl_file(os.path.abspath(filename), tree_file_abs_path, out_file, ctl_file,
                                       seq_type=seq_type, codon_frequency=codon_frequency, noisy=noisy, verbose=verbose,
                                       runmode=runmode, clock=clock, aminoacid_distance=aminoacid_distance, model=2,
                                       nssites=2, genetic_code=genetic_code, fix_kappa=fix_kappa,
                                       kappa=kappa, fix_omega=fix_omega_dict[model], omega=1, getSE=getSE,
                                       RateAncestor=RateAncestor, Mgene=0,
                                       small_difference=small_difference, clean_data=clean_data, method=method)
        self.parallel_execute(options_list, dir_list=dir_list)

        results_dict = OrderedDict()
        double_delta_dict = OrderedDict()
        raw_pvalues_dict = OrderedDict()
        raw_pvalues_list = []

        for basename in basename_dir_list:
            results_dict[basename] = OrderedDict()
            for model in model_list:
                output_file = "%s/%s/%s/%s.out" % (out_dir, basename, model, basename)
                codeml_report = CodeMLReport(output_file)
                results_dict[basename][model] = codeml_report.LnL

        for basename in basename_dir_list:
            for model in model_list:
                if results_dict[basename][model] is None:
                    print("LnL was not calculated for %s" % basename)
                    break
            else:
                doubled_delta = 2 * (results_dict[basename]["Model_A"] - results_dict[basename]["Model_A_null"])
                p_value = chisqprob(doubled_delta, 1)  # degrees of freedom = 1

                double_delta_dict[basename] = doubled_delta
                raw_pvalues_dict[basename] = p_value
                raw_pvalues_list.append(p_value)

        adjusted_pvalues_list = fdrcorrection0(raw_pvalues_list)
        i = 0
        with open(results_file, "w") as out_fd:
            out_fd.write("id\tmodel_a_null,LnL\tmodel_a,LnL\t2*delta\traw p-value\tadjusted p-value\n")
            for basename in basename_dir_list:
                for model in model_list:
                    if results_dict[basename][model] is None:
                        print("LnL was not calculated for %s" % basename)
                        i += 1
                        break
                else:
                    #doubled_delta = 2 * (results_dict[basename]["Model_A"] - results_dict[basename]["Model_A_null"])
                    #p_value = chisqprob(doubled_delta, 1) # degrees of freedom = 1
                    out_fd.write("%s\t%f\t%f\t%f\t%f\t%f\n" % (basename, results_dict[basename]["Model_A_null"],
                                                           results_dict[basename]["Model_A"], double_delta_dict[basename],
                                                           raw_pvalues_dict[basename], adjusted_pvalues_list[i]))
                    i += 1

