#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
from Tools.Abstract import Tool
from Routines import FileRoutines


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
                          genetic_code=0, fix_kappa=False, kappa=5, fix_omega=False, omega=0.2, getSE=0, RateAncestor=0,
                          small_difference=0.00001, clean_data=True, method=0):

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
        options += "runmode = %i\n\n" % runmode

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
        options += "omega = %f\n" % omega

        options += "getSE = %i\n" % getSE
        options += "RateAncestor = %i\n" % RateAncestor

        options += "Small_Diff = %f\n" % small_difference
        options += "cleandata = %i  * remove sites with ambiguity data (1:yes, 0:no)?\n" % (1 if clean_data else 0)
        options += "method = %i   * 0: simultaneous; 1: one branch at a time\n" % method

        with open(ctl_file, "w") as ctl_fd:
            ctl_fd.write(options)

    def run_codeml(self, ctl_file):

        options = " %s" % ctl_file

        self.execute(options=options)

    def parallel_codeml(self, in_dir, tree_file, out_dir, seq_type="codons", codon_frequency="F3X4", noisy=3,
                          verbose="concise", runmode=0, clock=0, aminoacid_distance=None, model=1, nssites=0,
                          genetic_code=0, fix_kappa=False, kappa=5, fix_omega=False, omega=0.2, getSE=0, RateAncestor=0,
                          small_difference=0.00001, clean_data=True, method=0):

        FileRoutines.save_mkdir(out_dir)
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
            FileRoutines.save_mkdir(filename_out_dir)
            self.generate_ctl_file(os.path.abspath(filename), tree_file_abs_path, out_file, ctl_file, seq_type=seq_type,
                                   codon_frequency=codon_frequency, noisy=noisy, verbose=verbose, runmode=runmode,
                                   clock=clock, aminoacid_distance=aminoacid_distance, model=model, nssites=nssites,
                                   genetic_code=genetic_code, fix_kappa=fix_kappa, kappa=kappa, fix_omega=fix_omega,
                                   omega=omega, getSE=getSE, RateAncestor=RateAncestor,
                                   small_difference=small_difference, clean_data=clean_data, method=method)
            self.parallel_execute(options_list, dir_list=dir_list)