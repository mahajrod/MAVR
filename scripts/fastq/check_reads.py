#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import FastQRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--forward_reads", action="store", dest="forward_reads", required=True,
                    help="Forward reads")
parser.add_argument("-r", "--reverse_reads", action="store", dest="reverse_reads", required=True,
                    help="Reverse reads")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

args = parser.parse_args()

with FastQRoutines.metaopen(args.forward_reads, "r") as in_fwd_fd, \
     FastQRoutines.metaopen(args.reverse_reads, "r") as in_rev_fd, \
     FastQRoutines.metaopen(f"{args.output_prefix}.clean_1.fastq", "w") as good_fwd_fd, \
     FastQRoutines.metaopen(f"{args.output_prefix}.clean_2.fastq", "w") as good_rev_fd, \
     FastQRoutines.metaopen(f"{args.output_prefix}.bad_1.fastq", "w") as bad_fwd_fd, \
     FastQRoutines.metaopen(f"{args.output_prefix}.bad_2.fastq", "w") as bad_rev_fd:

    for line in in_fwd_fd:
        fwd_name = line
        fwd_seq = in_fwd_fd.readline()
        fwd_sep = in_fwd_fd.readline()
        fwd_qual = in_fwd_fd.readline()

        rev_name = in_rev_fd.readline()
        rev_seq = in_rev_fd.readline()
        rev_sep = in_rev_fd.readline()
        rev_qual = in_rev_fd.readline()

        if (fwd_name[0] != "@") or (rev_name[0] != "@") or (fwd_sep[0] != "+") or (rev_sep[0] != "+"):
            for row in fwd_name, fwd_seq, fwd_sep, fwd_qual:
                bad_fwd_fd.write(row)
            for row in rev_name, rev_seq, rev_sep, rev_qual:
                bad_rev_fd.write(row)
        elif (len(fwd_seq) != len(fwd_qual)) or (len(rev_seq) != len(rev_qual)):
            for row in fwd_name, fwd_seq, fwd_sep, fwd_qual:
                bad_fwd_fd.write(row)
            for row in rev_name, rev_seq, rev_sep, rev_qual:
                bad_rev_fd.write(row)
        else:
            for row in fwd_name, fwd_seq, fwd_sep, fwd_qual:
                good_fwd_fd.write(row)
            for row in rev_name, rev_seq, rev_sep, rev_qual:
                good_rev_fd.write(row)