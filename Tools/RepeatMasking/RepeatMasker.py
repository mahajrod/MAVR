#!/usr/bin/env python
import os

from Tools.Abstract import Tool
from CustomCollections.GeneralCollections import IdSet
from Routines import FileRoutines


class RepeatMasker(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "repeatmasker", path=path, max_threads=max_threads)
        self.repeat_classes_used_for_gene_annotation = ["DNA",
                                                        "DNA\?",
                                                        "LINE",
                                                        "LINE\?",
                                                        "LTR",
                                                        "LTR\?",
                                                        "RC",
                                                        "RC\?",
                                                        "Retroposon",
                                                        "Retroposon\?",
                                                        "SINE",
                                                        "SINE\?",
                                                        "Helitron",
                                                        "Helitron\?"
                                                        ]

    @staticmethod
    def convert_rm_out_to_gff(input_file, output_file, annotated_repeat_classes_file, annotated_repeat_families_file):
        repeat_classes_set = IdSet()
        repeat_families_set = IdSet()
        with open(input_file, "r") as in_fd:
            for i in range(0, 3):
                in_fd.readline()

            with open(output_file, "w") as out_fd:
                for line in in_fd:
                    tmp = line.strip().split()
                    strand = "+" if tmp[8] == "+" else "-"
                    repeat_class_family = tmp[10].split("/")
                    if len(repeat_class_family) == 1:
                        repeat_class_family.append(".")
                    repeat_classes_set.add(repeat_class_family[0])
                    repeat_families_set.add("/".join(repeat_class_family))
                    parameters = "Class=%s;Family=%s;Matching_repeat=%s;SW_score=%s;Perc_div=%s;Perc_del=%s;Pers_ins=%s" \
                                 % (repeat_class_family[0], repeat_class_family[1],
                                    tmp[9], tmp[0], tmp[1], tmp[2], tmp[3])
                    out_fd.write("%s\tRepeatMasker\trepeat\t%s\t%s\t.\t%s\t.\t%s\n" % (tmp[4], tmp[5], tmp[6], strand, parameters))
        repeat_classes_set.write(annotated_repeat_classes_file)
        repeat_families_set.write(annotated_repeat_families_file)

    @staticmethod
    def extract_annotated_repeat_types_from_gff(gff_file, annotated_repeat_classes_file):
        sed_string = "sed -r 's/.*Class=(.*);Family.*/\1/' %s | sort | uniq > %s" % (gff_file,
                                                                                     annotated_repeat_classes_file)
        # awk variant of string
        # awk -F'\t' '{print $9}' repeatmasker.selected_repeat_classes.gff | awk -F';' '{print $1}' | awk -F'=' '{print $2}' | sort | uniq
        os.system(sed_string)

    def extract_repeats_used_for_gene_annotation(self, input_gff, output_gff):
        grep_pattern = "|".join(self.repeat_classes_used_for_gene_annotation)
        grep_string = "grep -P '%s'" % grep_pattern
        grep_string += " %s" % input_gff
        grep_string += " > %s" % output_gff
        self.execute(cmd=grep_string)

    def extract_repeats_from_database(self, output_file, species=None, clade=None, stat_mode=None):

        options = " -species %s" % species if species else ""
        options += " -clade %s" % clade if clade else ""
        options += " -stat" if stat_mode else ""
        options += " | awk 'NR > 1 {print $0}'" if not stat_mode else ""
        options += " > %s" % output_file

        self.execute(options=options, cmd="queryRepeatDatabase.pl")

    def mask(self, list_of_fasta_files, output_dir="./", soft_masking=True, engine="ncbi", slow_search=True,
             quick_search=False, rush_search=False, no_low_complexity=None, only_low_complexity=None,
             no_interspersed=None, only_interspersed=None, no_rna=None, only_alu=None, custom_library=None,
             species=None, html_output=False, ace_output=False, gff_output=False):

        if (slow_search and quick_search) or (rush_search and quick_search) or (slow_search and rush_search):
            raise ValueError("Both quick search(-q) and slow search(-s) options were set. Choose ONE!")

        if species and custom_library:
            tmp_repeat_file = "%s.repeats.tmp.fa" % species
            tmp_repeats_all_file = "all.repeats.tmp.fasta"
            self.extract_repeats_from_database(tmp_repeat_file, species=species)

            cmd = "cat %s %s > %s" % (tmp_repeat_file, custom_library, tmp_repeats_all_file)
            self.execute(cmd=cmd)

        options = " -pa %i" % self.threads
        options += " -e %s" % engine
        options += " -s" if slow_search else ""
        options += " -q" if quick_search else ""
        options += " -qq" if rush_search else ""
        options += " -nolow" if no_low_complexity else ""
        options += " -low" if only_low_complexity else ""
        options += " -noint" if no_interspersed else ""
        options += " -int" if only_interspersed else ""
        options += " -norna" if no_rna else ""
        options += " -alu" if only_alu else ""

        if species and custom_library:
            options += " -lib %s" % tmp_repeats_all_file
        elif custom_library:
            options += " -lib %s" % custom_library if custom_library else ""
        elif species:
            options += " -species %s" % species if species else ""

        options += " -dir %s" % output_dir
        options += " -html" if html_output else ""
        options += " -ace" if ace_output else ""
        options += " -gff" if gff_output else ""
        options += " -xsmall" if soft_masking else ""

        options += " " + (list_of_fasta_files if isinstance(list_of_fasta_files, str) else  " ".join(FileRoutines.make_list_of_path_to_files(list_of_fasta_files)))

        self.execute(options=options)

        """
    -source
        Includes for each annotation the HSP "evidence". Currently this
        option is only available with the "-html" output format listed
        below.

    -u  Creates an additional annotation file not processed by
        ProcessRepeats

    -xm Creates an additional output file in cross_match format (for
        parsing)



    -div [number]
        Masks only those repeats < x percent diverged from consensus seq


    -cutoff [number]
        Sets cutoff score for masking repeats when using -lib (default 225)



    Contamination options

    -is_only
        Only clips E coli insertion elements out of fasta and .qual files

    -is_clip
        Clips IS elements before analysis (default: IS only reported)

    -no_is
        Skips bacterial insertion element check

    Running options

    -gc [number]
        Use matrices calculated for 'number' percentage background GC level

    -gccalc
        RepeatMasker calculates the GC content even for batch files/small
        seqs

    -frag [number]
        Maximum sequence length masked without fragmenting (default 60000,
        300000 for DeCypher)

    -nocut
        Skips the steps in which repeats are excised

    -noisy
        Prints search engine progress report to screen (defaults to .stderr
        file)

    -nopost
        Do not postprocess the results of the run ( i.e. call ProcessRepeats
        ). NOTE: This options should only be used when ProcessRepeats will
        be run manually on the results.

    output options



    -a(lignments)
        Writes alignments in .align output file

    -inv
        Alignments are presented in the orientation of the repeat (with
        option -a)

    -lcambig
        Outputs ambiguous DNA transposon fragments using a lower case name.
        All other repeats are listed in upper case. Ambiguous fragments
        match multiple repeat elements and can only be called based on
        flanking repeat information.

    -small
        Returns complete .masked sequence in lower case


    -poly
        Reports simple repeats that may be polymorphic (in file.poly)


    -no_id
        Leaves out final column with unique ID for each element (was
        default)

    -e(xcln)
        Calculates repeat densities (in .tbl) excluding runs of >=20 N/Xs in
        the query
        """