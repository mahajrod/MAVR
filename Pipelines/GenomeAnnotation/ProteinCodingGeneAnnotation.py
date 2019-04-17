#!/usr/bin/env python
import os
from RouToolPa.Tools.BLAST import BLASTp
from RouToolPa.Collections.General import IdSet
from RouToolPa.Tools.Annotation import AUGUSTUS
from RouToolPa.Tools.Expression import Gffread
from RouToolPa.Tools.Alignment import STAR
from RouToolPa.Tools.HMMER import HMMER3
from Pipelines.Filtering import FilteringPipeline




class ProteinCodingGeneAnnotation(FilteringPipeline):

    def __init__(self):
        FilteringPipeline.__init__(self)

    def prepare_gene_annotation_directories(self, output_directory, protein_species_list=None,
                                            rnaseq_tissues_list=None):

        annotation_dir = "%s/annotation/" % output_directory
        masking_dir = "%s/masking/" % output_directory
        genome_dir = "%s/genome/" % output_directory

        protein_evidence_dir = "%s/protein_evidence/" % annotation_dir
        protein_evidence_fasta_dir = "%s/fasta/" % protein_evidence_dir
        protein_evidence_exonerate_dir = "%s/exonerate/" % protein_evidence_dir

        rnaseq_evidence_dir = "%s/rnaseq_evidence/" % annotation_dir
        rnaseq_read_dir = "%s/reads/" % rnaseq_evidence_dir
        rnaseq_alignment_dir = "%s/alignment/" % rnaseq_evidence_dir

        est_evidence_dir = "%s/est_evidence/" % annotation_dir

        augustus_dir = "%s/augustus/" % annotation_dir
        augustus_hints_dir = "%s/hints/" % augustus_dir
        augustus_splited_hints_dir = "%s/splited/" % augustus_hints_dir
        augustus_raw_dir = "%s/raw/" % augustus_dir
        augustus_db_or_hint_support_dir = "%s/supported_db_or_hint/" % augustus_dir
        augustus_db_and_hint_support_dir = "%s/supported_db_and_hint/" % augustus_dir

        augustus_pfam_alignments = "%s/pfam/" % augustus_dir
        augustus_swissprot_alignments = "%s/swissprot/" % augustus_dir

        for directory in (genome_dir, masking_dir,
                          annotation_dir, protein_evidence_dir, protein_evidence_fasta_dir,
                          protein_evidence_exonerate_dir, rnaseq_evidence_dir, rnaseq_read_dir,
                          rnaseq_alignment_dir, est_evidence_dir, augustus_dir, augustus_hints_dir,
                          augustus_splited_hints_dir, augustus_raw_dir, augustus_db_or_hint_support_dir,
                          augustus_db_and_hint_support_dir, augustus_pfam_alignments,
                          augustus_swissprot_alignments):
            self.safe_mkdir(directory)

        for species in protein_species_list:
            self.safe_mkdir("%s/%s/" % (protein_evidence_dir, species))
        for tissue in rnaseq_tissues_list:
            self.safe_mkdir("%s/%s/" % (rnaseq_evidence_dir, tissue))
            #for sample in sample_list:
            #    FileRoutines.save_mkdir("%s/%s" % (directory, sample))

        return (genome_dir, masking_dir,
                annotation_dir, protein_evidence_dir, protein_evidence_fasta_dir,
                protein_evidence_exonerate_dir, rnaseq_evidence_dir, rnaseq_read_dir,
                rnaseq_alignment_dir, est_evidence_dir, augustus_dir, augustus_hints_dir,
                augustus_splited_hints_dir, augustus_raw_dir, augustus_db_or_hint_support_dir,
                augustus_db_and_hint_support_dir, augustus_pfam_alignments,
                augustus_swissprot_alignments)

    def softmask_genome(self, nonmasked_genome, repeatmasker_out_file):
        pass

    def prepare_rnaseq_hints(self, rnaseq_read_dir, rnaseq_alignment_dir, genome_dir, hints_dir, genome_fasta=None,
                             samples=None, annotation_gtf=None, sjdboverhang=None,
                             genomeSAindexNbases=None, genomeChrBinNbits=None, genome_size=None,
                             feature_from_gtf_to_use_as_exon=None, exon_tag_to_use_as_transcript_id=None,
                             exon_tag_to_use_as_gene_id=None, length_of_sequences_flanking_junction=None,
                             junction_tab_file_list=None, three_prime_trim=None, five_prime_trim=None,
                             adapter_seq_for_three_prime_clip=None, max_mismatch_percent_for_adapter_trimming=None,
                             three_prime_trim_after_adapter_clip=None, output_type="BAM", sort_bam=True,
                             max_memory_for_bam_sorting=8000000000, include_unmapped_reads_in_bam=True,
                             output_unmapped_reads=True, two_pass_mode=True, max_intron_length=None,
                             min_reads_supporting_junction_hint=1, source="RNASEQ", priority=100, threads=1,
                             STAR_path=""):
        STAR.threads = threads
        STAR.path = STAR_path

        samples = samples if samples else self.get_sample_list(rnaseq_read_dir)

        STAR.align_samples(rnaseq_read_dir, rnaseq_alignment_dir, genome_dir, genome_fasta=genome_fasta,
                           samples=samples,
                           annotation_gtf=annotation_gtf,  sjdboverhang=sjdboverhang,
                           genomeSAindexNbases=genomeSAindexNbases, genomeChrBinNbits=genomeChrBinNbits,
                           genome_size=genome_size, feature_from_gtf_to_use_as_exon=feature_from_gtf_to_use_as_exon,
                           exon_tag_to_use_as_transcript_id=exon_tag_to_use_as_transcript_id,
                           exon_tag_to_use_as_gene_id=exon_tag_to_use_as_gene_id,
                           length_of_sequences_flanking_junction=length_of_sequences_flanking_junction,
                           junction_tab_file_list=junction_tab_file_list, three_prime_trim=three_prime_trim,
                           five_prime_trim=five_prime_trim,
                           adapter_seq_for_three_prime_clip=adapter_seq_for_three_prime_clip,
                           max_mismatch_percent_for_adapter_trimming=max_mismatch_percent_for_adapter_trimming,
                           three_prime_trim_after_adapter_clip=three_prime_trim_after_adapter_clip,
                           output_type=output_type, sort_bam=sort_bam,
                           max_memory_for_bam_sorting=max_memory_for_bam_sorting,
                           include_unmapped_reads_in_bam=include_unmapped_reads_in_bam,
                           output_unmapped_reads=output_unmapped_reads, two_pass_mode=two_pass_mode,
                           max_intron_length=max_intron_length)

        for sample in samples:
            sample_junction_file = "%s/%s/SJ.out.tab " % (rnaseq_alignment_dir, sample)
            sample_hint_gff = "%s/rnaseq.hints.%s.gff" % (hints_dir, sample)
            AUGUSTUS.convert_star_junctions_to_intron_hints(sample_junction_file, sample_hint_gff,
                                                            min_supporting_reads=min_reads_supporting_junction_hint,
                                                            source=source, priority=priority)

    def prepare_protein_hints(self):
        pass

    def prepare_est_hints(self):
        pass

    def prepare_hints(self):
        pass

    def predict_genes(self, output_prefix, annotation_species_prefix, genome_fasta, augustus_species, output_directory="./", augustus_strand=None,
                      augustus_gene_model=None, augustus_config_dir=None, augustus_use_softmasking=None,
                      augustus_other_options="", augustus_hintsfile=None, augustus_extrinsicCfgFile=None,
                      augustus_predict_UTR=None, augustus_min_intron_len=None, threads=1, augustus_dir="",
                      hmmer_dir="", blast_dir="",
                      stop_codons_list=("TGA", "TAA", "TAG"), genetic_code_table=1):

        draft_file_prefix = "%s/raw/%s" % (output_directory, output_prefix)

        augustus_splited_input_dir = "%s/splited_input/" % output_directory
        augustus_splited_output_dir = "%s/splited_output_dir" % output_directory

        output_raw_gff = "%s.raw.gff" % draft_file_prefix
        output_gff = "%s.renamed.gff" % draft_file_prefix
        augustus_pep = "%s.pep" % draft_file_prefix

        AUGUSTUS.path = augustus_dir
        AUGUSTUS.threads = threads
        HMMER3.path = hmmer_dir
        HMMER3.threads = threads
        BLASTp.path = blast_dir
        BLASTp.threads = threads

        print("Annotating genes...")
        AUGUSTUS.parallel_predict(augustus_species, genome_fasta, output_raw_gff, strand=augustus_strand,
                                  gene_model=augustus_gene_model, output_gff3=True,
                                  other_options=augustus_other_options, config_dir=augustus_config_dir,
                                  use_softmasking=augustus_use_softmasking, hints_file=augustus_hintsfile,
                                  split_dir=augustus_splited_input_dir, splited_output_dir=augustus_splited_output_dir,
                                  extrinsicCfgFile=augustus_extrinsicCfgFile, predict_UTR=augustus_predict_UTR,
                                  combine_output_to_single_file=True, min_intron_len=augustus_min_intron_len)

        #replace_augustus_ids(augustus_gff, output_prefix, species_prefix=None, number_of_digits_in_id=8):

        AUGUSTUS.replace_augustus_ids(output_raw_gff, draft_file_prefix, species_prefix=annotation_species_prefix,
                                      number_of_digits_in_id=8)
        #extract_transcript_sequences(self, input_gff_file, genomic_fasta_file, output_prefix, coding_only=False)
        gffread_file_prefix = "%s.gffread" % draft_file_prefix
        gffread_transcripts_file, gffread_cds_file, gffread_pep_file = Gffread.extract_transcript_sequences(output_gff,
                                                                                                            genome_fasta,
                                                                                                            gffread_file_prefix)
        gffread_trimmed_cds = ".".join(gffread_cds_file.split(".")[:-1]) + ".trimmed.cds"
        gffread_trimmed_pep = ".".join(gffread_pep_file.split(".")[:-1]) + ".trimmed.pep"
        self.trim_cds_and_remove_terminal_stop_codons(gffread_cds_file, gffread_trimmed_cds,
                                                      stop_codons_list=stop_codons_list) # using default stop_codons(from universal genetic_code)/ Note that this will affect mtDNA proteins
        inframe_stop_codons_file_prefix = "%s.inframe_stop_codon" % draft_file_prefix
        self.translate_sequences_from_file(gffread_trimmed_cds, gffread_trimmed_pep,
                                           format="fasta", id_expression=None,
                                           genetic_code_table=genetic_code_table, translate_to_stop=False,
                                           prefix_of_file_inframe_stop_codons_seqsin=inframe_stop_codons_file_prefix) # Universal code !!!

        AUGUSTUS.extract_gene_ids_from_output(output_gff, all_annotated_genes_ids)
        AUGUSTUS.extract_CDS_annotations_from_output(output_gff, CDS_gff)

        print("Extracting peptides...")

        AUGUSTUS.extract_proteins_from_output(output_gff, output_pep, id_prefix="", evidence_stats_file=output_evidence_stats,
                                              supported_by_hints_file=output_supported_stats)

        self.compare_sequences_from_files(output_pep, "%s.trimmed.pep" % args.output, "comparison_of_peptides",
                                                      format="fasta", verbose=True)

        os.system("awk -F'\\t' 'NR==1 {}; NR > 1 {print $2}' %s > %s" % (output_supported_stats, output_supported_stats_ids))

        print("Annotating domains(Pfam database)...")

        HMMER3.parallel_hmmscan(args.pfam_db, output_pep, output_hmmscan, num_of_seqs_per_scan=None,
                                split_dir="splited_hmmscan_fasta/",
                                splited_output_dir="splited_hmmscan_output_dir",
                                tblout_outfile=None, domtblout_outfile=output_domtblout, pfamtblout_outfile=None,
                                splited_tblout_dir=None, splited_domtblout_dir="hmmscan_domtblout/")
        HMMER3.extract_dom_ids_hits_from_domtblout(output_domtblout, output_pfam_annotated_dom_ids)
        hits_dict = HMMER3.extract_dom_names_hits_from_domtblout(output_domtblout, output_pfam_annotated_dom_names)
        supported_ids = IdSet(hits_dict.keys())
        supported_ids.write(output_pfam_supported_transcripts_ids)
        remove_transcript_ids_str = "sed -re 's/\.t[0123456789]+//' %s | sort -k 1 | uniq > %s" % (output_pfam_supported_transcripts_ids,
                                                                                                   output_pfam_supported_genes_ids)
        os.system(remove_transcript_ids_str)


        print("Annotating peptides(Swissprot database)...")

        BLASTp.parallel_blastp(output_pep, args.swissprot_db, evalue=0.0000001, output_format=6,
                               outfile=output_swissprot_blastp_hits, split_dir="splited_blastp_fasta",
                               splited_output_dir="splited_blastp_output_dir")
        hits_dict = BLASTp.extract_hits_from_tbl_output(output_swissprot_blastp_hits, output_swissprot_blastp_hits_names)
        supported_ids = IdSet(hits_dict.keys())
        supported_ids.write(output_swissprot_supported_transcripts_ids)

        remove_transcript_ids_str = "sed -re 's/\.t[0123456789]+//' %s | sort -k 1 | uniq > %s" % (output_swissprot_supported_transcripts_ids,
                                                                                                   output_swissprot_supported_genes_ids)
        os.system(remove_transcript_ids_str)




        """

output_evidence_stats = "%s.transcript.evidence" % args.output
output_evidence_stats_longest_pep = "%s.transcript.evidence.longest_pep" % args.output
output_supported_stats = "%s.transcript.supported" % args.output
output_supported_stats_ids = "%s.transcript.supported.ids" % args.output
output_supported_stats_longest_pep = "%s.transcript.supported.longest_pep" % args.output
output_hmmscan = "%s.hmmscan.hits" % args.output
output_domtblout = "%s.domtblout" % args.output
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


for directory in ("splited_blastp_fasta", "splited_blastp_output_dir"):
    shutil.rmtree(directory)

gene_ids_black_list = [genes_masked_ids] if args.masking else []
gene_ids_white_list = []


gene_ids_white_list = [output_pfam_supported_genes_ids, output_swissprot_supported_genes_ids]
HMMER3.intersect_ids_from_files(output_swissprot_supported_transcripts_ids, output_pfam_supported_transcripts_ids,
                                output_swissprot_pfam_supported_transcripts_ids, mode="combine")
HMMER3.intersect_ids_from_files(output_swissprot_pfam_supported_transcripts_ids, output_supported_stats_ids,
                                output_swissprot_pfam_or_hints_supported_transcripts_ids, mode="combine")
HMMER3.intersect_ids_from_files(output_swissprot_pfam_supported_transcripts_ids, output_supported_stats_ids,
                                output_swissprot_pfam_and_hints_supported_transcripts_ids, mode="common")

SequenceRoutines.extract_sequence_by_ids(output_pep, output_swissprot_pfam_or_hints_supported_transcripts_ids,
                                         output_swissprot_pfam_or_hints_supported_transcripts_pep)
SequenceRoutines.extract_sequence_by_ids(output_pep, output_swissprot_pfam_and_hints_supported_transcripts_ids,
                                         output_swissprot_pfam_and_hints_supported_transcripts_pep)

AUGUSTUS.extract_evidence_by_ids(output_evidence_stats, output_swissprot_pfam_or_hints_supported_transcripts_ids,
                                 output_swissprot_pfam_or_hints_supported_transcripts_evidence)
AUGUSTUS.extract_evidence_by_ids(output_evidence_stats, output_swissprot_pfam_and_hints_supported_transcripts_ids,
                                 output_swissprot_pfam_and_hints_supported_transcripts_evidence)
AUGUSTUS.extract_longest_isoforms(output_swissprot_pfam_or_hints_supported_transcripts_evidence,
                                  output_swissprot_pfam_or_hints_supported_transcripts_longest_pep_evidence)
AUGUSTUS.extract_longest_isoforms(output_swissprot_pfam_and_hints_supported_transcripts_evidence,
                                  output_swissprot_pfam_and_hints_supported_transcripts_longest_pep_evidence)

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

HMMER3.intersect_ids_from_files([all_annotated_genes_ids], gene_ids_black_list, genes_not_masked_ids, mode="only_a")
HMMER3.intersect_ids_from_files(gene_ids_white_list, gene_ids_black_list, final_genes_ids, mode="only_a")

final_ids = IdSet()
final_ids.read(final_genes_ids)

AnnotationsRoutines.extract_annotation_from_gff(output_gff, final_ids, ["gene"], final_gff)
AUGUSTUS.extract_CDS_annotations_from_output(final_gff, final_CDS_gff)

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

if args.pfam_db and args.swissprot_db:
    db_or_hints_dir = "supported_by_db_or_hints/"
    db_and_hints_dir = "supported_by_db_and_hints/"
    for directory in db_and_hints_dir, db_or_hints_dir:
        FileRoutines.save_mkdir(directory)

    os.system("mv %s.supported.transcripts.swissprot_or_pfam_or_hints* %s" % (args.output, db_or_hints_dir))
    os.system("mv %s.supported.transcripts.swissprot_or_pfam_and_hints* %s" % (args.output, db_and_hints_dir))

"""

