#!/usr/bin/env python3
__author__ = 'mahajrod'
import argparse
from subprocess import PIPE, Popen
from RouToolPa.Tools.Bedtools import Intersect
from RouToolPa.Parsers.FastQC import FastQCReport



parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fastqc", action="store", dest="fastqc",
                    help="fastqc report file")
parser.add_argument("-o", "--output", action="store", dest="out",
                    help="output bed file")
parser.add_argument("-b", "--blastdb", action="store", dest="blastdb",
                    help="blast database")
parser.add_argument("-t", "--blast_task", action="store", dest="blast_task", default="blastn-short",
                    help="blast mode")
parser.add_argument("-l", "--logfile", action="store", dest="logfile", default="log.log",
                    help="log")
parser.add_argument("-g", "--gff", action="store", dest="gfffile",
                    help="gff with annotations")
parser.add_argument("-c", "--cutoff", action="store", dest="cutoff", type=float, default=0.01,
                    help="e-value cut-off")
parser.add_argument("-i", "--intersection_file", action="store", dest="inter_file",
                    help="reulting intersection file", default="intersection.bed")
parser.add_argument("-u", "--unfound", action="store", dest="unfound_file",
                    help="unfound sequences", default="unfound.t")
args = parser.parse_args()


fastqc_report = FastQCReport(fastqc_report_file=args.fastqc)
overrepresented_sequences = fastqc_report.overrepresented_sequences()

index = 1
bed_string_list = []
found_sequences = set([])
unfound_sequences = set([])
print("Handling sequences...\nOnly hits with e_value <= %f are shown" % args.cutoff)
with open(args.logfile, "w") as log_fd:
    log_fd.write("blast_fastqc.py logfile\n\nHandling sequences...\n")
    for sequence in overrepresented_sequences:
        blastn_output = Popen(["echo \"%s\" | blastn -task %s -num_threads 8 -db %s -outfmt \"10 sacc sstart send sstrand length  qstart qend qlen evalue\" "
                               % (sequence, args.blast_task, args.blastdb)], shell=True, stdout=PIPE)
        print("\n" + sequence)
        log_fd.write("Sequence:\n\t%s\nHits:\n" % sequence)
        blastn_output_list = [line for line in blastn_output.stdout]
        found_flag = False
        for hit in blastn_output_list:
            hit_name = "Hit_%i" % index
            hit_string = hit.decode('utf-8').strip()
            log_fd.write("\t%s\n" % hit_string)
            hit_string_list = hit_string.split(",")
            if hit_string_list[3] == "minus":
                hit_string_list[1], hit_string_list[2] = hit_string_list[2], hit_string_list[1]
            e_value = float(hit_string_list[8])
            if e_value <= args.cutoff:
                #print(hit_string_list)
                bed_string_list.append("\t".join(hit_string_list[:3] +
                                                 [hit_name, hit_string_list[8], hit_string_list[3], sequence])
                                       + "\n")

                print(hit_string)
                found_flag = True
                index += 1
                continue

        if found_flag:
            found_sequences.add(sequence)
        else:
            unfound_sequences.add(sequence)
    log_fd.write("\nSummary:\n\tTotaly\t%i\n\tFound\t%i\n" % (len(overrepresented_sequences), len(found_sequences)))
    for sequence in found_sequences:
        log_fd.write("\t\t%s\n" % sequence)
    log_fd.write("\tUnfound\t%i\n" % (len(unfound_sequences)))
    for sequence in unfound_sequences:
        log_fd.write("\t\t%s\n" % sequence)

with open(args.out, "w") as bed_fd:
    for bed_string in bed_string_list:
        bed_fd.write(bed_string)
with open(args.unfound_file, "w") as uf_fd:
    for uf_sequence in unfound_sequences:
        uf_fd.write(uf_sequence)

if args.gfffile:
    Intersect.intersect(args.gfffile, args.out, args.inter_file, method="-u")
