#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
import argparse
from RouToolPa.Tools.BLAST import BLASTp
from RouToolPa.Tools.HMMER import HMMER3
from RouToolPa.Tools.Bedtools import Intersect
from RouToolPa.Tools.Annotation import AUGUSTUS
from RouToolPa.Tools.Expression import Gffread
from RouToolPa.Routines import AnnotationsRoutines, MatplotlibRoutines, SequenceRoutines, FileRoutines
from RouToolPa.Collections.General import IdSet




parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", required=True,
                    help="Input fasta file with sequences")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Prefix of output files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads to use")
parser.add_argument("-x", "--species_prefix", action="store", dest="species_prefix", required=True,
                    help="Species-related prefix to use in IDs")
parser.add_argument("-s", "--species", action="store", dest="species", required=True,
                    help="Species to use as model")
parser.add_argument("-r", "--strand", action="store", dest="strand", default="both",
                    help="Strand to consider. Possible variants: both, forward, backward."
                         "Default: both")
parser.add_argument("-g", "--gene_model", action="store", dest="gene_model",
                    help="Gene model to use. Possible variants:"
                         "partial      : allow prediction of incomplete genes at the sequence boundaries (default)"
                         "intronless   : only predict single-exon genes like in prokaryotes and some eukaryotes"
                         "complete     : only predict complete genes"
                         "atleastone   : predict at least one complete gene"
                         "exactlyone   : predict exactly one complete gene")
parser.add_argument("-e", "--other_options", action="store", dest="other_options",
                    help="Other augustus options")
parser.add_argument("-c", "--augustus_config_dir", action="store", dest="config_dir",
                    help="Augustus config dir")
parser.add_argument("-p", "--pfam_hmm3", action="store", dest="pfam_db",
                    help="Pfam database in hmm3 format")
parser.add_argument("-w", "--swissprot_blast_db", action="store", dest="swissprot_db",
                    help="Blast database of swissprot")
parser.add_argument("-m", "--masking", action="store", dest="masking",
                    help="Gff of bed file with masking of repeats")
parser.add_argument("--softmasking", action="store_true", dest="softmasking",
                    help="Use softmasking from genome")
parser.add_argument("--hintsfile", action="store", dest="hintsfile",
                    help="File with hints")
parser.add_argument("--extrinsicCfgFile", action="store", dest="extrinsicCfgFile",
                    help="Config file with scoring for hints")
parser.add_argument("-u", "--predict_UTR", action="store_true", dest="predict_UTR",
                    help="Predict UTR. works not for all species")
parser.add_argument("-a", "--augustus_dir", action="store", dest="augustus_dir", default="",
                    help="Directory with augustus binary")

args = parser.parse_args()

output_raw_gff = "%s.gff" % args.output
output_gff = "%s.renamed.gff" % args.output
output_pep = "%s.pep" % args.output
output_evidence_stats = "%s.transcript.evidence" % args.output
output_evidence_stats_longest_pep = "%s.transcript.evidence.longest_pep" % args.output
output_supported_stats = "%s.transcript.supported" % args.output
output_supported_stats_ids = "%s.transcript.supported.ids" % args.output
output_supported_stats_longest_pep = "%s.transcript.supported.longest_pep" % args.output
output_prefix_hmmscan = "%s.hmmscan" % args.output
output_domtblout = "%s.hmmscan.domtblout" % args.output
output_pfam_annotated_dom_ids = "%s.pfam.dom_ids" % args.output
#output_pfam_supported_ids = "%s.supported.pfam.ids" % args.output
output_pfam_supported_transcripts_ids = "%s.supported.transcripts.pfam.ids" % args.output
output_pfam_supported_genes_ids = "%s.supported.genes.pfam.ids" % args.output

output_pfam_annotated_dom_names = "%s.pfam.dom_names" % args.output

output_swissprot_blastp_hits = "%s.swissprot.hits" % args.output
#output_swissprot_supported_ids = "%s.supported.swissprot.ids" % args.output
output_swissprot_supported_transcripts_ids = "%s.supported.transcripts.swissprot.ids" % args.output
output_swissprot_supported_genes_ids = "%s.supported.genes.swissprot.ids" % args.output
output_swissprot_blastp_hits_names = "%s.swissprot.hits.names" % args.output

output_swissprot_pfam_supported_transcripts_ids = "%s.supported.transcripts.swissprot_or_pfam.ids" % args.output
output_swissprot_pfam_or_hints_supported_transcripts_ids = "%s.supported.transcripts.swissprot_or_pfam_or_hints.ids" % args.output
output_swissprot_pfam_or_hints_supported_transcripts_inframe_stop_ids = "%s.supported.transcripts.swissprot_or_pfam_or_hints.inframe_stop.ids" % args.output
output_swissprot_pfam_and_hints_supported_transcripts_ids = "%s.supported.transcripts.swissprot_or_pfam_and_hints.ids" % args.output
output_swissprot_pfam_and_hints_supported_transcripts_inframe_stop_ids = "%s.supported.transcripts.swissprot_or_pfam_and_hints.inframe_stop.ids" % args.output
output_swissprot_pfam_or_hints_supported_transcripts_evidence = "%s.supported.transcripts.swissprot_or_pfam_or_hints.evidence" % args.output
output_swissprot_pfam_and_hints_supported_transcripts_evidence = "%s.supported.transcripts.swissprot_or_pfam_and_hints.evidence" % args.output
output_swissprot_pfam_or_hints_supported_transcripts_pep = "%s.supported.transcripts.swissprot_or_pfam_or_hints.pep" % args.output
output_swissprot_pfam_and_hints_supported_transcripts_pep = "%s.supported.transcripts.swissprot_or_pfam_and_hints.pep" % args.output

output_swissprot_pfam_or_hints_supported_transcripts_longest_pep_evidence = "%s.supported.transcripts.swissprot_or_pfam_or_hints.longest_pep.evidence" % args.output
output_swissprot_pfam_and_hints_supported_transcripts_longest_pep_evidence = "%s.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.evidence" % args.output
output_swissprot_pfam_or_hints_supported_transcripts_longest_pep = "%s.supported.transcripts.swissprot_or_pfam_or_hints.longest_pep.pep" % args.output
output_swissprot_pfam_and_hints_supported_transcripts_longest_pep = "%s.supported.transcripts.swissprot_or_pfam_and_hints.longest_pep.pep" % args.output

prefix_of_file_inframe_stop_codons_seqs = "cds_with_in_frame_stop_codons"
cds_with_inframe_stop_codons_ids = "%s.ids" % prefix_of_file_inframe_stop_codons_seqs

CDS_gff = "%s.CDS.gff" % args.output
CDS_masked_gff = "%s.CDS.masked.gff" % args.output
all_annotated_genes_ids = "%s.genes.all.ids" % args.output
genes_masked_ids = "%s.genes.masked.ids" % args.output
genes_not_masked_ids = "%s.genes.not.masked.ids" % args.output
final_genes_ids = "%s.genes.final.ids" % args.output

final_gff = "%s.final.gff" % args.output
final_CDS_gff = "%s.final.CDS.gff" % args.output

AUGUSTUS.path = args.augustus_dir
AUGUSTUS.threads = args.threads

print("Annotating genes...")

AUGUSTUS.parallel_predict(args.species, args.input, output_raw_gff, strand=args.strand, gene_model=args.gene_model,
                          output_gff3=True, other_options=args.other_options, config_dir=args.config_dir,
                          use_softmasking=args.softmasking, hints_file=args.hintsfile,
                          extrinsicCfgFile=args.extrinsicCfgFile, predict_UTR=args.predict_UTR,
                          parsing_mode="parse")

AUGUSTUS.replace_augustus_ids(output_raw_gff, args.output, species_prefix=args.species_prefix,
                              number_of_digits_in_id=8)

Gffread.extract_transcript_sequences(output_gff, args.input, args.output)

SequenceRoutines.trim_cds_and_remove_terminal_stop_codons("%s.cds" % args.output, "%s.trimmed.cds" % args.output,
                                                          stop_codons_list=("TGA", "TAA", "TAG")) # using default stop_codons(from universal genetic_code)/ Note that this will affect mtDNA proteins

SequenceRoutines.translate_sequences_from_file("%s.trimmed.cds" % args.output, "%s.trimmed.pep" % args.output,
                                               format="fasta", id_expression=None,
                                               genetic_code_table=1, translate_to_stop=False,
                                               prefix_of_file_inframe_stop_codons_seqs=prefix_of_file_inframe_stop_codons_seqs) # Universal code !!!

AUGUSTUS.extract_gene_ids_from_output(output_gff, all_annotated_genes_ids)
AUGUSTUS.extract_CDS_annotations_from_output(output_gff, CDS_gff)
if args.masking:
    print("Intersecting annotations with repeats...")
    Intersect.intersect(CDS_gff, args.masking, CDS_masked_gff, method="-u")
    sed_string = "sed 's/.*=//;s/\.t.*//' %s | sort | uniq > %s" % (CDS_masked_gff, genes_masked_ids)
    os.system(sed_string)

print("Extracting peptides...")

AUGUSTUS.extract_proteins_from_output(output_gff, output_pep, id_prefix="", evidence_stats_file=output_evidence_stats,
                                      supported_by_hints_file=output_supported_stats)

SequenceRoutines.compare_sequences_from_files(output_pep, "%s.trimmed.pep" % args.output, "comparison_of_peptides",
                                              format="fasta", verbose=True)

os.system("awk -F'\\t' 'NR==1 {}; NR > 1 {print $2}' %s > %s" % (output_supported_stats, output_supported_stats_ids))

if args.pfam_db:
    print("Annotating domains(Pfam database)...")
    HMMER3.threads = args.threads
    HMMER3.parallel_hmmscan(args.pfam_db, output_pep, output_prefix_hmmscan, "./", #output_hmmscan,
                            num_of_seqs_per_scan=None) # TODO CHECK!!!!!!!!!!!!!!!!!!!!!!!!!!
    HMMER3.extract_dom_ids_hits_from_domtblout(output_domtblout, output_pfam_annotated_dom_ids)
    hits_dict = HMMER3.extract_dom_names_hits_from_domtblout(output_domtblout, output_pfam_annotated_dom_names)
    supported_ids = IdSet(hits_dict.keys())
    supported_ids.write(output_pfam_supported_transcripts_ids)
    remove_transcript_ids_str = "sed -re 's/\.t[0123456789]+//' %s | sort -k 1 | uniq > %s" % (output_pfam_supported_transcripts_ids,
                                                                                               output_pfam_supported_genes_ids)
    os.system(remove_transcript_ids_str)

if args.swissprot_db:
    print("Annotating peptides(Swissprot database)...")
    BLASTp.threads = args.threads
    BLASTp.parallel_blastp(output_pep, args.swissprot_db, evalue=0.0000001, output_format=6,
                           outfile=output_swissprot_blastp_hits, split_dir="splited_blastp_fasta",
                           splited_output_dir="splited_blastp_output_dir")
    hits_dict = BLASTp.extract_hits_from_tbl_output(output_swissprot_blastp_hits, output_swissprot_blastp_hits_names)
    supported_ids = IdSet(hits_dict.keys())
    supported_ids.write(output_swissprot_supported_transcripts_ids)

    remove_transcript_ids_str = "sed -re 's/\.t[0123456789]+//' %s | sort -k 1 | uniq > %s" % (output_swissprot_supported_transcripts_ids,
                                                                                               output_swissprot_supported_genes_ids)
    os.system(remove_transcript_ids_str)

    for directory in ("splited_blastp_fasta", "splited_blastp_output_dir"):
        shutil.rmtree(directory)

gene_ids_black_list = [genes_masked_ids] if args.masking else []
gene_ids_white_list = []

print("Combining database results...")
if args.pfam_db and args.swissprot_db:

    gene_ids_white_list = [output_pfam_supported_genes_ids, output_swissprot_supported_genes_ids]
    HMMER3.intersect_ids_from_files(output_swissprot_supported_transcripts_ids, output_pfam_supported_transcripts_ids,
                                    output_swissprot_pfam_supported_transcripts_ids, mode="combine")
    HMMER3.intersect_ids_from_files(output_swissprot_pfam_supported_transcripts_ids, output_supported_stats_ids,
                                    output_swissprot_pfam_or_hints_supported_transcripts_ids, mode="combine")
    HMMER3.intersect_ids_from_files(output_swissprot_pfam_supported_transcripts_ids, output_supported_stats_ids,
                                    output_swissprot_pfam_and_hints_supported_transcripts_ids, mode="common")

    print("Extracting sequences...")
    SequenceRoutines.extract_sequence_by_ids(output_pep, output_swissprot_pfam_or_hints_supported_transcripts_ids,
                                             output_swissprot_pfam_or_hints_supported_transcripts_pep)
    SequenceRoutines.extract_sequence_by_ids(output_pep, output_swissprot_pfam_and_hints_supported_transcripts_ids,
                                             output_swissprot_pfam_and_hints_supported_transcripts_pep)

    print("Extracting evidence...")
    AUGUSTUS.extract_evidence_by_ids(output_evidence_stats, output_swissprot_pfam_or_hints_supported_transcripts_ids,
                                     output_swissprot_pfam_or_hints_supported_transcripts_evidence)
    AUGUSTUS.extract_evidence_by_ids(output_evidence_stats, output_swissprot_pfam_and_hints_supported_transcripts_ids,
                                     output_swissprot_pfam_and_hints_supported_transcripts_evidence)
    print("Extracting longest isoforms...")
    AUGUSTUS.extract_longest_isoforms(output_swissprot_pfam_or_hints_supported_transcripts_evidence,
                                      output_swissprot_pfam_or_hints_supported_transcripts_longest_pep_evidence)
    AUGUSTUS.extract_longest_isoforms(output_swissprot_pfam_and_hints_supported_transcripts_evidence,
                                      output_swissprot_pfam_and_hints_supported_transcripts_longest_pep_evidence)
    print("Extracting evidence...")
    SequenceRoutines.extract_sequence_by_ids(output_pep, "%s.ids" % output_swissprot_pfam_or_hints_supported_transcripts_longest_pep_evidence,
                                             output_swissprot_pfam_or_hints_supported_transcripts_longest_pep)
    SequenceRoutines.extract_sequence_by_ids(output_pep, "%s.ids" % output_swissprot_pfam_and_hints_supported_transcripts_longest_pep_evidence,
                                             output_swissprot_pfam_and_hints_supported_transcripts_longest_pep)

    for id_file in output_swissprot_pfam_or_hints_supported_transcripts_ids, \
                   output_swissprot_pfam_and_hints_supported_transcripts_ids, \
                   "%s.ids" % output_swissprot_pfam_or_hints_supported_transcripts_longest_pep_evidence, \
                   "%s.ids" % output_swissprot_pfam_and_hints_supported_transcripts_longest_pep_evidence:
        out_pref = id_file[:-4]
        out_gff = "%s.gff" % out_pref
        AnnotationsRoutines.extract_transcripts_by_ids(output_gff, id_file, out_gff)
        for suffix in ".trimmed.cds", ".transcript":
            SequenceRoutines.extract_sequence_by_ids("%s%s" % (args.output, suffix),
                                                     id_file,
                                                     "%s%s" % (out_pref, suffix))

    HMMER3.intersect_ids_from_files(output_swissprot_pfam_or_hints_supported_transcripts_ids,
                                    cds_with_inframe_stop_codons_ids,
                                    output_swissprot_pfam_or_hints_supported_transcripts_inframe_stop_ids,
                                    mode="common")

    HMMER3.intersect_ids_from_files(output_swissprot_pfam_and_hints_supported_transcripts_ids,
                                    cds_with_inframe_stop_codons_ids,
                                    output_swissprot_pfam_and_hints_supported_transcripts_inframe_stop_ids,
                                    mode="common")
elif args.pfam_db:
    gene_ids_white_list = [output_pfam_supported_genes_ids]
elif args.swissprot_db:
    gene_ids_white_list = [output_swissprot_supported_genes_ids]
else:
    gene_ids_white_list = [all_annotated_genes_ids]

HMMER3.intersect_ids_from_files([all_annotated_genes_ids], gene_ids_black_list, genes_not_masked_ids, mode="only_a")
HMMER3.intersect_ids_from_files(gene_ids_white_list, gene_ids_black_list, final_genes_ids, mode="only_a")
#"""
final_ids = IdSet()
final_ids.read(final_genes_ids)


AnnotationsRoutines.extract_annotation_from_gff(output_gff, final_ids, ["gene"], final_gff)
print("Extracting CDS...")
AUGUSTUS.extract_CDS_annotations_from_output(final_gff, final_CDS_gff)

print("Drawing histograms...")

for stat_file in output_evidence_stats, output_supported_stats, \
                 output_swissprot_pfam_or_hints_supported_transcripts_longest_pep_evidence, \
                 output_swissprot_pfam_and_hints_supported_transcripts_longest_pep_evidence, \
                 output_swissprot_pfam_or_hints_supported_transcripts_evidence, \
                 output_swissprot_pfam_and_hints_supported_transcripts_evidence:

    MatplotlibRoutines.percent_histogram_from_file(stat_file, stat_file, data_type=None,
                                                   column_list=(2,),
                                                   comments="#", n_bins=20,
                                                   title="Transcript support by hints",
                                                   extensions=("png", "svg"),
                                                   legend_location="upper center",
                                                   stats_as_legend=True)
print("Creating final directories...")
if args.pfam_db and args.swissprot_db:
    db_or_hints_dir = "supported_by_db_or_hints/"
    db_and_hints_dir = "supported_by_db_and_hints/"
    for directory in db_and_hints_dir, db_or_hints_dir:
        FileRoutines.safe_mkdir(directory)

    os.system("mv %s.supported.transcripts.swissprot_or_pfam_or_hints* %s" % (args.output, db_or_hints_dir))
    os.system("mv %s.supported.transcripts.swissprot_or_pfam_and_hints* %s" % (args.output, db_and_hints_dir))
