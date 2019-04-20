#/usr/bin/env python
import os
import collections
from RouToolPa.Routines import FileRoutines





from collections import OrderedDict

from Bio import SeqIO

from RouToolPa.Tools.AssemblyTools import spades
#from RouToolPa.Tools.FilterTools import trim_galore
#from RouToolPa.Tools.AlignmentTools import Bowtie2, BWA
from RouToolPa.Tools.Alignment import Bowtie2, BWA, Novoalign, TMAP

from Pipelines.Abstract import Pipeline

from RouToolPa.Tools.Bedtools import GenomeCov
from RouToolPa.Tools.Samtools import SamtoolsV1
from RouToolPa.Tools.Picard import MarkDuplicates
from RouToolPa.Tools.GATK import *
from RouToolPa.Tools.Picard import CreateSequenceDictionary

from RouToolPa.Collections.General import SynDict, IdSet

# TODO: refactor, remove absolete functions and so on


class SNPCallPipeline(Pipeline):

    def __init__(self, max_threads=1, max_memory=10, GATK_dir="", GATK_jar="GenomeAnalysisTK.jar",
                 Picard_dir=""):
        Pipeline.__init__(self, max_threads=max_threads, max_memory=max_memory)

        self.GATK_dir = GATK_dir
        self.GATK_jar = GATK_jar
        self.Picard_dir = Picard_dir

    def filter_reference(self, reference, repeatmasking_gff_list, output_prefix, reference_len_file=None,
                         annotation_gff=None, max_masked_fraction=0.8, white_scaffold_list=(), black_scaffold_list=(),
                         max_length=None):
        scaffold_with_annotation_file = "%s.scaffolds_with_annotations.ids" % output_prefix
        sorted_combined_repeatmasking_gff = "%s.sorted_combined_repeatmasking.gff" % output_prefix
        reference_len_filename = "%s.reference.len" % output_prefix if reference_len_file is None else reference_len_file
        repeatmasking_coverage_file = "%s.repeatmasking.coverage" % output_prefix
        filtering_log_file = "%s.filtering.log" % output_prefix

        filtered_scaffolds_file = "%s.filtered.fasta" % output_prefix
        filtered_out_scaffolds_file = "%s.filtered_out.fasta" % output_prefix
        filtered_out_scaffolds_id_file = "%s.filtered_out.ids" % output_prefix

        scaffold_with_annotation_set = self.get_scaffold_ids_from_gff(annotation_gff,
                                                                      out_file=scaffold_with_annotation_file)

        print("Sorting GFFs with masking...")
        sorting_string = "cat %s | sort -k1,1 -k4,4n -k5,5n > %s" % (repeatmasking_gff_list if isinstance(repeatmasking_gff_list, str) else" ".join(repeatmasking_gff_list),
                                                                     sorted_combined_repeatmasking_gff)

        self.execute(options=sorting_string, cmd="")

        print("Parsing reference...")

        reference_dict = self.parse_seq_file(reference, mode="parse")

        if reference_len_file is None:
            length_dict = self.get_lengths(reference_dict, out_file=reference_len_filename)
        else:
            length_dict = SynDict(filename=reference_len_filename)

        print("Calculating fraction of masked regions...")

        GenomeCov.get_coverage_for_gff(sorted_combined_repeatmasking_gff, reference_len_filename,
                                       output=repeatmasking_coverage_file)

        print("Filtering...")
        low_zero_coverage_fraction_dict = SynDict(filename=repeatmasking_coverage_file, key_index=0, value_index=4,
                                                  include_line_expression=lambda l: l.split("\t")[1] == "0",
                                                  expression=float,
                                                  include_value_expression=lambda v: v < (1.0 - max_masked_fraction))
        #print low_zero_coverage_fraction_dict
        scaffold_to_remove = IdSet()

        with open(filtering_log_file, "w") as log_fd:
            log_fd.write("#Scaffold\tStatus\tLength\tDescription\n")
            for scaffold in reference_dict:
                if scaffold in black_scaffold_list:
                    scaffold_to_remove.add(scaffold)
                    log_fd.write("%s\tRemoved\t%i\tBlackList\n" % (scaffold, length_dict[scaffold]))
                    continue

                if scaffold in low_zero_coverage_fraction_dict:
                    if scaffold in white_scaffold_list:
                        log_fd.write("%s\tRetained\t%i\tWhiteList,LowNonMaskedPercentage:%f\n" % (scaffold, length_dict[scaffold], low_zero_coverage_fraction_dict[scaffold]))
                        continue
                    if scaffold in scaffold_with_annotation_set:
                        log_fd.write("%s\tRetained\t%i\tWithAnnotations,LonNonMaskedPercentage:%f\n" % (scaffold, length_dict[scaffold], low_zero_coverage_fraction_dict[scaffold]))
                        continue
                    if not(max_length is None):
                        if length_dict[scaffold] > max_length:
                            log_fd.write("%s\tRetained\t%i\tLong,LowNonMaskedPercentage:%f\n" % (scaffold, length_dict[scaffold], low_zero_coverage_fraction_dict[scaffold]))
                        continue

                    scaffold_to_remove.add(scaffold)
                    log_fd.write("%s\tRemoved\t%i\tLowNonMaskedPercentage:%f\n" % (scaffold, length_dict[scaffold], low_zero_coverage_fraction_dict[scaffold]))
                    continue
                log_fd.write("%s\tRetained\t%i\tOK\n" % (scaffold, length_dict[scaffold]))

        scaffold_to_remove.write(filename=filtered_out_scaffolds_id_file)

        SeqIO.write(self.record_by_id_generator(reference_dict, scaffold_to_remove),
                    filtered_out_scaffolds_file, format="fasta")
        SeqIO.write(self.record_by_id_generator(reference_dict, IdSet(reference_dict.keys()) - scaffold_to_remove),
                    filtered_scaffolds_file, format="fasta")
        print("Total scaffolds\t%i\nRemoved\t%i\n" % (len(reference_dict), len(scaffold_to_remove)))

    def prepare_dirs(self, sample_list, outdir="./", include_alignment_dir=False):
        dir_dict = {
                    "SNPcall/": OrderedDict()
                    }

        if include_alignment_dir:
            dir_dict["alignment/"] = OrderedDict()

        for subdir in dir_dict:
            for sample in sample_list:
                dir_dict[subdir][sample] = OrderedDict()
                
        self.recursive_mkdir(dir_dict, out_dir=outdir)

    def call_variants(self, sample_dir, reference, merged_prefix, sample_list=None, outdir="./", suffix=None, input="alignment",
                      input_filetype="bam", threads=None, mark_duplicates=False, known_variants_vcf=None,
                      genotyping_mode="DISCOVERY", output_mode="EMIT_VARIANTS_ONLY",
                      stand_call_conf=30, skip_base_score_recalibration=False,
                      iteration_number=3, SNP_QD=2.0, SNP_FS=30.0, SNP_MQ=40.0, SNP_MappingQualityRankSum=-12.5,
                      SNP_ReadPosRankSum=-8.0, indel_QD=2.0, indel_ReadPosRankSum=-20.0, indel_FS=200.0,
                      SNP_filter_name="ambiguous_snp", indel_filter_name="ambiguous_indel",
                      analyze_covariates=True, include_region_id_file=None, exclude_region_id_file=None):

        SamtoolsV1.check_for_fasta_index(reference)

        CreateSequenceDictionary.jar_path = self.Picard_dir
        CreateSequenceDictionary.check_for_fasta_dict(reference)

        for tool in VariantFiltration, \
                    MarkDuplicates, \
                    RealignerTargetCreator, \
                    IndelRealigner, \
                    BaseRecalibrator, \
                    PrintReads, \
                    HaplotypeCaller, \
                    SelectVariants, \
                    GenotypeGVCFs,\
                    CombineVariants:

            tool.threads = threads if threads else self.threads
            tool.max_memory = "%ig" % self.max_memory
            tool.jar_path = self.GATK_dir

        samples = self.get_sample_list(sample_dir, sample_list=sample_list)

        self.prepare_dirs(samples, outdir=outdir, include_alignment_dir=input == "reads")

        if input == "reads":
            pass
        elif input == "alignment":
            alignment_filename_prefix_template = "%%s%s" % suffix

        known_sites = known_variants_vcf

        iterations = 1 if skip_base_score_recalibration else iteration_number # do only one(zero) iteration if skip base recalibration

        for sample in samples:
            if mark_duplicates:
                """
                java -Xmx100g -jar ~/tools/picard-tools-2.5.0/picard.jar MarkDuplicates \
                     I=${bam} \
                     O=${bam%bam}rmdup.bam \
                      M=${bam}.mark_dup_metrics.txt

                java -jar ~/tools/picard-tools-2.5.0/picard.jar BuildBamIndex \
                        INPUT=${bam%bam}rmdup.bam

                """

                pass
                # sample_alignment_prefix =
            else:
                sample_alignment_prefix = "%s/%s/%s" % (sample_dir, sample, alignment_filename_prefix_template % sample)

            sample_alignment = "%s.%s" % (sample_alignment_prefix, input_filetype)
            sample_intervals_for_realignment = "%s.forIndelRealigner.intervals" % sample_alignment_prefix
            sample_realigned_bam = "%s.realigned.bam" % sample_alignment_prefix

            """
            RealignerTargetCreator.create(reference, sample_alignment,
                                          output=sample_intervals_for_realignment,
                                          known_indels_vcf=None,
                                          max_interval_size=None,
                                          min_reads_cov=None,
                                          mismatch_fraction=None,
                                          window_size=None,
                                          default_base_qualities=None)

            IndelRealigner.realign(reference,
                                   sample_alignment,
                                   sample_realigned_bam,
                                   target_intervals=sample_intervals_for_realignment,
                                   known_indels_vcf=None, model=None, lod_threshold=None,
                                   entropy_threshold=None, max_cons=None,
                                   max_size_for_movement=None, max_pos_move=None, max_reads_for_cons=None,
                                   max_reads_for_realignment=None, max_reads_in_memory=None, no_original_tags=False,
                                   nway_out=False, default_base_qualities=None)
            """
        for iteration_index in range(0, iterations):
            gvcf_list = []

            sample_recall_table = "%s.recall_data.iteration%i.grp" % (sample_alignment_prefix, iteration_index)
            sample_postrecall_table = "%s.postrecall_data.iteration%i.grp" % (sample_alignment_prefix, iteration_index)
            sample_recall_plots = "%s.recall.iteration%i.pdf" % (sample_alignment_prefix, iteration_index)
            sample_recall_csv = "%s.recall.iteration%i.csv" % (sample_alignment_prefix, iteration_index)

            sample_recalled_reads_bam = "%s.recal_reads.iteration%i.bam" % (sample_alignment_prefix, iteration_index)

            merged_vcf_prefix = "%s/SNPcall/%s.iteration%i" % (outdir, merged_prefix, iteration_index)
            merged_raw_vcf_prefix = "%s.raw" % merged_vcf_prefix
            merged_raw_vcf = "%s.vcf" % merged_raw_vcf_prefix

            merged_raw_snp_vcf = "%s.raw.snp.vcf" % merged_vcf_prefix
            merged_with_filters_snp_vcf = "%s.with_filters.snp.vcf" % merged_vcf_prefix
            merged_filtered_snp_vcf = "%s.filtered.snp.vcf" % merged_vcf_prefix

            merged_raw_indel_vcf = "%s.raw.indel.vcf" % merged_vcf_prefix
            merged_with_filters_indel_vcf = "%s.with_filters.indel.vcf" % merged_vcf_prefix
            merged_filtered_indel_vcf = "%s.filtered.indel.vcf" % merged_vcf_prefix

            merged_filtered_combined_vcf = "%s.filtered.combined.vcf" % merged_vcf_prefix

            for sample in samples:
                vcf_prefix = "%s/SNPcall/%s/%s.iteration%i" % (outdir, sample, sample, iteration_index)
                gvcf = "%s.g.vcf" % vcf_prefix
                #raw_snp_vcf = "%s.raw.snp.gvcf" % vcf_prefix
                #raw_indel_vcf = "%s.raw.indel.gvcf" % vcf_prefix

                sample_alignment_prefix = "%s/%s/%s" % (sample_dir, sample, alignment_filename_prefix_template % sample)
                sample_realigned_bam = "%s.realigned.bam" % sample_alignment_prefix

                gvcf_list.append(gvcf)

                if ((not skip_base_score_recalibration) and known_sites is not None) or (iteration_index > 0):

                    BaseRecalibrator.get_recalibration_table(reference,
                                                             sample_realigned_bam,
                                                             sample_recall_table,
                                                             known_sites,
                                                             include_region_id_file=include_region_id_file,
                                                             exclude_region_id_file=exclude_region_id_file)

                    BaseRecalibrator.get_recalibration_table(reference,
                                                             sample_realigned_bam,
                                                             sample_postrecall_table,
                                                             known_sites,
                                                             BQSR=sample_recall_table,
                                                             include_region_id_file=include_region_id_file,
                                                             exclude_region_id_file=exclude_region_id_file)
                    if analyze_covariates:
                        AnalyzeCovariates.plot_two_recall_table(reference, sample_recall_table, sample_postrecall_table,
                                                                sample_recall_plots, csv_out=sample_recall_csv)

                    PrintReads.get_recalled_reads(reference,
                                                  sample_realigned_bam,
                                                  sample_recall_table,
                                                  sample_recalled_reads_bam)

                    #HaplotypeCaller.call(reference, sample_realigned, raw_vcf, genotyping_mode=genotyping_mode,
                    #                     output_mode=output_mode, stand_emit_conf=stand_emit_conf, stand_call_conf=stand_call_conf)
                    """
                    HaplotypeCaller.gvcf_call(reference, sample_recalled_reads_bam, gvcf, genotyping_mode=genotyping_mode,
                                              output_mode=output_mode, stand_call_conf=stand_call_conf,
                                              include_region_id_file=include_region_id_file,
                                              exclude_region_id_file=exclude_region_id_file)
                    """
                else:

                    HaplotypeCaller.gvcf_call(reference, sample_realigned_bam, gvcf, genotyping_mode=genotyping_mode,
                                              output_mode=output_mode, stand_call_conf=stand_call_conf,
                                              include_region_id_file=include_region_id_file,
                                              exclude_region_id_file=exclude_region_id_file)

            GenotypeGVCFs.genotype(reference, gvcf_list, merged_raw_vcf_prefix)

            self.hardfilter_variants(reference, merged_raw_vcf, merged_vcf_prefix, SNP_QD=SNP_QD, SNP_FS=SNP_FS,
                                     SNP_MQ=SNP_MQ, SNP_MappingQualityRankSum=SNP_MappingQualityRankSum,
                                     SNP_ReadPosRankSum=SNP_ReadPosRankSum, indel_QD=indel_QD,
                                     indel_ReadPosRankSum=indel_ReadPosRankSum, indel_FS=indel_FS,
                                     SNP_filter_name=SNP_filter_name, indel_filter_name=indel_filter_name,
                                     threads=threads)

            """
            SelectVariants.select_variants(reference, merged_raw_vcf, merged_raw_snp_vcf, vartype="SNP", varfilter=None)
            SelectVariants.select_variants(reference, merged_raw_vcf, merged_raw_indel_vcf, vartype="INDEL", varfilter=None)

            VariantFiltration.filter_bad_SNP(reference, merged_raw_snp_vcf, merged_with_filters_snp_vcf,
                                             filter_name=SNP_filter_name,
                                             QD=SNP_QD, FS=SNP_FS, MQ=SNP_MQ,
                                             MappingQualityRankSum=SNP_MappingQualityRankSum,
                                             ReadPosRankSum=SNP_ReadPosRankSum)
            VariantFiltration.filter_bad_indel(reference, merged_raw_indel_vcf, merged_with_filters_indel_vcf,
                                               filter_name=indel_filter_name, QD=indel_QD,
                                               ReadPosRankSum=indel_ReadPosRankSum, FS=indel_FS)

            SelectVariants.remove_entries_with_filters(reference, merged_with_filters_snp_vcf, merged_filtered_snp_vcf)
            SelectVariants.remove_entries_with_filters(reference, merged_with_filters_indel_vcf, merged_filtered_indel_vcf)

            CombineVariants.combine_from_same_source(reference,
                                                     [merged_filtered_snp_vcf, merged_filtered_indel_vcf],
                                                     merged_filtered_combined_vcf)
            """
            known_sites = merged_filtered_combined_vcf

    def hardfilter_variants(self, reference, raw_vcf, output_prefix, SNP_QD=2.0, SNP_FS=30.0, SNP_MQ=40.0,
                            SNP_MappingQualityRankSum=-12.5,
                            SNP_ReadPosRankSum=-8.0, indel_QD=2.0, indel_ReadPosRankSum=-20.0, indel_FS=200.0,
                            SNP_filter_name="ambiguous_snp", indel_filter_name="ambiguous_indel", threads=None):

        raw_snp_vcf = "%s.raw.snp.vcf" % output_prefix
        raw_indel_vcf = "%s.raw.indel.vcf" % output_prefix
        with_filters_snp_vcf = "%s.with_filters.snp.vcf" % output_prefix
        with_filters_indel_vcf = "%s.with_filters.indel.vcf" % output_prefix
        filtered_snp_vcf = "%s.filtered.snp.vcf" % output_prefix
        filtered_indel_vcf = "%s.filtered.indel.vcf" % output_prefix
        filtered_combined_vcf = "%s.filtered.combined.vcf" % output_prefix

        for tool in VariantFiltration, \
                    SelectVariants, \
                    CombineVariants:

            tool.threads = threads if threads else self.threads
            tool.max_memory = "%ig" % self.max_memory
            tool.jar_path = self.GATK_dir

        SelectVariants.select_variants(reference, raw_vcf, raw_snp_vcf, vartype="SNP", varfilter=None)
        SelectVariants.select_variants(reference, raw_vcf, raw_indel_vcf, vartype="INDEL", varfilter=None)

        VariantFiltration.filter_bad_SNP(reference, raw_snp_vcf, with_filters_snp_vcf,
                                         filter_name=SNP_filter_name,
                                         QD=SNP_QD, FS=SNP_FS, MQ=SNP_MQ,
                                         MappingQualityRankSum=SNP_MappingQualityRankSum,
                                         ReadPosRankSum=SNP_ReadPosRankSum)
        VariantFiltration.filter_bad_indel(reference, raw_indel_vcf, with_filters_indel_vcf,
                                           filter_name=indel_filter_name, QD=indel_QD,
                                           ReadPosRankSum=indel_ReadPosRankSum, FS=indel_FS)

        SelectVariants.remove_entries_with_filters(reference, with_filters_snp_vcf, filtered_snp_vcf)
        SelectVariants.remove_entries_with_filters(reference, with_filters_indel_vcf, filtered_indel_vcf)

        CombineVariants.combine_from_same_source(reference,
                                                 [filtered_snp_vcf, filtered_indel_vcf],
                                                 filtered_combined_vcf)


def get_chromosomes_bed(reference, reference_index, mitochondrial_region_name="mt",
                        chrom_out_file="chromosomes.bed", mito_out_file="mt.bed", reference_filetype="fasta"):
    if isinstance(reference, collections.MutableSequence):
        ref = reference
    else:
        ref = [reference]
    record_dict = SeqIO.index_db(reference_index, ref, reference_filetype)
    lengthes_dict = {}

    for record_id in record_dict:
        if record_id == mitochondrial_region_name:
            with open(mito_out_file, "w") as mt_fd:
                mt_fd.write(record_id + "\t1\t" + str(len(record_dict[record_id])) + "\n")
            continue
        lengthes_dict[record_id] = len(record_dict[record_id])

    with open(chrom_out_file, "w") as in_fd:
        for record_id in sorted(list(lengthes_dict.keys())):
            in_fd.write(record_id + "\t1\t" + str(lengthes_dict[record_id]) + "\n")


def alignment_sorting_and_filtering(sample_name,
                                    alignment_file,
                                    chromosomes_bed_file=None,
                                    mitochondrial_bed_file=None,
                                    input_alignment_type="sam",
                                    remove_SA=False,
                                    remove_nonPA=False,
                                    qualimap_stat=False):
    unaligned = 4
    supplementary_alignment = 2048
    non_primary_alignment = 256
    #-F 4 - skip UNaligned reads
    if input_alignment_type == "sam":
        input_type = "-S"
    elif input_alignment_type == "bam":
        input_type = ""

    os.system("samtools view %s -b %s | samtools sort - %s_trimmed_sorted" % (input_type, alignment_file, sample_name))
    os.system("samtools rmdup %s_trimmed_sorted.bam %s_trimmed_sorted_rm_pcr.bam" % (sample_name, sample_name))

    os.system("samtools view -b -f 4  %s_trimmed_sorted_rm_pcr.bam  > %s_trimmed_sorted_rm_pcr_unaligned.bam"
              % (sample_name, sample_name))

    extract_chromosomes = ""
    if chromosomes_bed_file:
        extract_chromosomes = "-L %s" % chromosomes_bed_file

    # -F 256 remove non primary alignment
    # -F 2048 remove supplementary alignment
    filter_parameter = unaligned
    if remove_SA:
        filter_parameter += supplementary_alignment
    if remove_nonPA:
        filter_parameter += non_primary_alignment

    os.system("samtools view -b -F %i %s %s_trimmed_sorted_rm_pcr.bam  > %s_trimmed_sorted_rm_pcr_chrom.bam"
              % (filter_parameter, extract_chromosomes, sample_name, sample_name))
    os.system("samtools index %s_trimmed_sorted_rm_pcr_chrom.bam" % sample_name)
    if qualimap_stat:
        os.system("qualimap bamqc -bam %s_trimmed_sorted_rm_pcr_chrom.bam " % sample_name)

    if mitochondrial_bed_file:
        os.system("samtools view -b -F %i -L %s %s_trimmed_sorted_rm_pcr.bam  > %s_trimmed_sorted_rm_pcr_mt.bam"
                  % (filter_parameter, mitochondrial_bed_file, sample_name, sample_name))
        os.system("samtools index %s_trimmed_sorted_rm_pcr_mt.bam" % sample_name)
        if qualimap_stat:
            os.system("qualimap bamqc -bam %s_trimmed_sorted_rm_pcr_mt.bam " % sample_name)

    os.system("rm -rf %s_trimmed.sam %s_trimmed_sorted.bam %s_trimmed_sorted_rm_pcr.bam"
              % (sample_name, sample_name, sample_name))


def get_alignment(index,
                  sample_name,
                  forward_reads,
                  chromosomes_bed_file=None,
                  mitochondrial_bed_file=None,
                  reverse_reads=None,
                  max_threads=5,
                  quality_score="phred33",
                  aligner="bowtie2",
                  remove_SA=False,
                  remove_nonPA=False,
                  find_discordant_alignments=True,
                  find_separated_alignments=True,
                  reads_filetype="fastq",
                  reference=None,
                  qualimap_stat=False):

    print("Handling %s sample..." % sample_name)
    output_file = "%s_trimmed.sam" % sample_name
    if aligner == "bowtie2":
        #bowtie2 = Bowtie2()
        Bowtie2.threads = max_threads
        Bowtie2.align(index, forward_reads, reverse_reads,
                      quality_score=quality_score,
                      alignment_mode="very-sensitive", output_file=output_file,
                      find_discordant_alignments=find_discordant_alignments,
                      find_separated_alignments=find_separated_alignments)
    elif aligner == "bwa-mem":
        #bwa = BWA()
        BWA.threads = max_threads
        BWA.mem(index, forward_reads, reverse_reads, output_file=output_file)
    elif aligner == "novoalign":
        # TODO: fix format choose
        #print(reads_filetype)
        #in_format = "FA" if reads_filetype == "fasta" else "ILM1.8"
        Novoalign.max_threads = max_threads
        Novoalign.align([forward_reads, reverse_reads] if reverse_reads else forward_reads,
                        index,
                        in_fmt="FA" if reads_filetype == "fasta" else "STDFQ",
                        output_file=output_file
                        )
    elif aligner == "tmap":
        TMAP.threads = max_threads
        TMAP.map4(reference, forward_reads, output=output_file)

    """
    if reverse_reads:
        os.system("bowtie2 --very-sensitive --%s -p %i -x %s -1 %s -2 %s > %s_trimmed.sam"
                  % (quality_score, max_threads, bowtie2_index, forward_reads, reverse_reads, sample_name))
    else:
        os.system("bowtie2 --very-sensitive --%s -p %i -x %s -U %s > %s_trimmed.sam"
                  % (quality_score, max_threads, bowtie2_index, forward_reads, sample_name))
    """
    alignment_sorting_and_filtering(sample_name, output_file, chromosomes_bed_file, mitochondrial_bed_file,
                                    remove_SA=remove_SA, remove_nonPA=remove_nonPA, qualimap_stat=qualimap_stat)


def filter_data(sample_name,
                min_length,
                forward_reads,
                forward_trim,
                reverse_reads=None,
                reverse_trim=None,
                quality_score="phred33",
                adapter="AGATCGGAAGAGC",
                quality_treshold=20,
                output_folder="trimmed",
                filter_tool="trim_galore"):
    print("Handling %s sample..." % sample_name)
    os.system("mkdir -p %s" % output_folder)

    if filter_tool == "trim_galore":
        trim_galore(min_length,
                    forward_reads,
                    forward_trim,
                    reverse_reads=reverse_reads,
                    reverse_trim=reverse_trim,
                    quality_score=quality_score,
                    adapter=adapter,
                    quality_treshold=quality_treshold,
                    output_folder=output_folder)

def assembly_reads(forward_reads,
                  reverse_reads=None,
                  max_threads=4,
                  platform="illumina",
                  output_dir="spades",
                  kmer_length_list=[]):
    #print("Handling %s sample..." % sample_name)
    spades(forward_reads, reverse_reads=reverse_reads,
           max_threads=max_threads, platform=platform,
           gzip_corrected_reads=False, output_dir=output_dir,
           kmer_length_list=kmer_length_list)


def get_coverage_thresholds(coverage_dist_file, one_side_base_threshold=0.025, minimum_threshold=5):
    #coverage_dist_file - is file like qualimap coverage_histogram.txt derived from alignment statistics

    fd = open(coverage_dist_file, "r")
    fd.readline()
    fd.readline()
    coverage = []
    frequency = []
    for line in fd:
        striped = line.strip()
        if striped == "":
            break
        #print (line)
        striped = striped.split("\t")
        coverage.append(float(striped[0]))
        frequency.append(float(striped[1]))
    fd.close()
    number_of_basses = sum(frequency)
    low_tr = int(one_side_base_threshold * float(number_of_basses))
    high_tr = int((1.00 - one_side_base_threshold) * float(number_of_basses))
    i = 0
    freq = 0
    while freq < low_tr:
        freq += frequency[i]
        i += 1
    min_coverage = coverage[i]
    while freq < high_tr:
        freq += frequency[i]
        i += 1
    max_coverage = coverage[i]
    return int(max(min_coverage, minimum_threshold)), int(max_coverage)


def snp_call(alignment,
             sample_name,
             reference_file,
             min_coverage,
             max_coverage,
             alignment_quality=40,
             snp_quality=100):

    os.system("samtools mpileup  -q %i -ugf %s %s | bcftools view -cvgN - > %s_raw.vcf"
              % (alignment_quality, reference_file, alignment,  sample_name))
    os.system("cat '%s_raw.vcf' | vcfutils.pl varFilter -D %i -d %i > %s.vcf"
              % (sample_name, max_coverage, min_coverage, sample_name))

    os.system("vcftools --vcf %s.vcf --out %s_filtered --remove-indels --recode --recode-INFO-all --minQ %i"
              % (sample_name, sample_name, snp_quality))


def snp_call_GATK(alignment,
                 sample_name,
                 reference_file,
                 known_sites_vcf,
                 stand_emit_conf=40,
                 stand_call_conf=100,
                 QD=2.0,
                 FS=60.0,
                 MQ=40.0,
                 HaplotypeScore=13.0,
                 MappingQualityRankSum=-12.5,
                 ReadPosRankSum=-8.0,
                 GATK_dir="",
                 num_of_threads=5,
                 skip_base_recalibration=False):
    #default filter expression
    #"QD < 2.0 || FS > 60.0 || MQ < 40.0 || HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 || ReadPosRankSum < -8.0"
    gatk_dir = FileRoutines.check_path(GATK_dir)
    intermediate_alignment = alignment
    if not skip_base_recalibration:
        intermediate_alignment = alignment + "_recal_reads.bam"
        #Analyze patterns of covariation in the sequence dataset
        os.system("java -jar %sGenomeAnalysisTK.jar -nct %i  -T BaseRecalibrator -R %s -I %s -knownSites %s -o %s_recal_data.table"
                  % (gatk_dir, num_of_threads, reference_file, alignment, known_sites_vcf, sample_name))
        #Do a second pass to analyze covariation remaining after recalibration
        os.system("java -jar %sGenomeAnalysisTK.jar -nct %i  -T BaseRecalibrator -R %s -I %s -knownSites %s  -BQSR %s_recal_data.table -o %s_post_recal_data.table"
                  % (gatk_dir, num_of_threads, reference_file, alignment, known_sites_vcf, sample_name, sample_name))

        #Generate before/after plots
        #os.system("java -jar %sGenomeAnalysisTK.jar -T AnalyzeCovariates -R %s -before %s_recal_data.table -after %s_post_recal_data.table -plots %s_recalibration_plots.pdf"
        #          % (gatk_dir, reference_file, sample_name, sample_name, sample_name))

        #Apply the recalibration to your sequence data
        os.system("java -jar %sGenomeAnalysisTK.jar -nct %i -T PrintReads -R %s -I %s -BQSR %s_recal_data.table -o %s"
                  % (gatk_dir, num_of_threads, reference_file, alignment, sample_name, intermediate_alignment))
    print("\nSNP call...\n")
    #SNP call
    os.system(" java -jar %sGenomeAnalysisTK.jar -nt %i -l INFO -R %s -T UnifiedGenotyper -I %s -stand_call_conf %i -stand_emit_conf %i  -o %s_GATK_raw.vcf --output_mode EMIT_VARIANTS_ONLY"
              % (gatk_dir, num_of_threads, reference_file, intermediate_alignment, stand_call_conf, stand_emit_conf, sample_name))
    #extract SNP
    os.system("java -jar %sGenomeAnalysisTK.jar -T SelectVariants -R %s -V %s_GATK_raw.vcf -selectType SNP -o %s_GATK_raw_no_indel.vcf"
              % (gatk_dir, reference_file, sample_name,  sample_name))
    #extract indels
    os.system("java -jar %sGenomeAnalysisTK.jar -T SelectVariants -R %s -V %s_GATK_raw.vcf -selectType INDEL -o %s_GATK_raw_only_indel.vcf"
              % (gatk_dir, reference_file, sample_name,  sample_name))

    #filtering
    print("\nFiltering SNP...\n")
    os.system("java -jar %sGenomeAnalysisTK.jar -T VariantFiltration -R %s -V %s_GATK_raw_no_indel.vcf --filterExpression 'QD < %f || FS > %f || MQ < %f || HaplotypeScore > %f || MappingQualityRankSum < %f || ReadPosRankSum < %f' --filterName 'ambigious_snp' -o %s_GATK_filtered_snps.vcf "
             % (gatk_dir, reference_file, sample_name, QD, FS, MQ, HaplotypeScore, MappingQualityRankSum, ReadPosRankSum, sample_name))
    #os.system("vcftools --vcf %s_GATK_filtered_snps.vcf --remove-filtered-all --out %s_GATK_best_snps.vcf --recode --recode-INFO-all"
    #          % (sample_name, sample_name ))

    """
    os.system("java -jar %sGenomeAnalysisTK.jar -nt %i -T HaplotypeCaller -R %s -I recal_reads.bam --genotyping_mode DISCOVERY --min_base_quality_score %i -stand_emit_conf %i -stand_call_conf %i -o %s"
              % (gatk_dir, num_of_threads, reference_file, min_base_quality_score, stand_emit_conf, stand_call_conf, raw_vcf_outfile))
    """


def snp_call_pipeline(bowtie2_index,
                      sample_name,
                      min_length,
                      reference_file,
                      reference_index,
                      right_reads,
                      right_trim,
                      left_reads=None,
                      left_trim=None,
                      skip_correction=False,
                      coverage_one_side_base_threshold=0.025,
                      coverage_minimum_threshold=5,
                      alignment_quality=40,
                      snp_quality=40,
                      max_threads=5,
                      chromosomes_bed_file="chromosomes.bed",
                      mitochondrial_bed_file="mt.bed",
                      mitochondrial_region_name="mt"):

    get_chromosomes_bed(reference, reference_index, mitochondrial_region_name=mitochondrial_region_name,
                        chrom_out_file=chromosomes_bed_file, mito_out_file=mitochondrial_bed_file,
                        reference_filetype="fasta")

    get_alignment(bowtie2_index, sample_name, min_length, right_reads,
                  right_trim, chromosomes_bed_file,
                  mitochondrial_bed_file, left_reads=left_reads, left_trim=left_trim,
                  skip_correction=skip_correction, max_threads=max_threads)

    for region in ["mt", "chrom"]:
        min_coverage, max_coverage = \
            get_coverage_thresholds("%s_trimmed_sorted_rm_pcr_%s_stats/raw_data/coverage_histogram.txt"
                                    % (sample_name, region),
                                    one_side_base_threshold=coverage_one_side_base_threshold,
                                    minimum_threshold=coverage_minimum_threshold)
        snp_call("%s_trimmed_sorted_rm_pcr_%s.bam" % (sample_name, region),
                 sample_name + "_" + region,
                 reference_file,
                 min_coverage,
                 max_coverage,
                 alignment_quality=alignment_quality,
                 snp_quality=snp_quality)


if __name__ == "__main__":
    samples_list =      [
                        "210-AID_Can1", #
                        "210-AID_Can2",    #
                        "210-Can1",  #
                        "210-Can2",   #
                        "210-FOA1",   #
                        "210-FOA2",
                        "210-Glu-Can2",    #1
                        "210-Glu-FOA2",    #
                        "210-Glu-FOA3",    #
                        "210-L1",   #
                        "210-L2",    #
                        "210-L3",   #
                        "210-L4",   #
                        "210-L5",  #
                        "210-L6",   #
                        "Sample_1",
                        "Sample_2",
                        "Sample_3",
                        "Sample_4",
                        "Sample_5",
                        "Sample_6",
                        "Sample_7",
                        "Sample_8",
                        "Sample_9",
                        "Sample_10",
                        "Sample_11",
                        "Sample_12",
                        "Sample_13",
                        "Sample_14",
                        "Sample_15",
                        "Sample_16",
                        "Sample_17",
                        "Sample_18",
                        "Sample_19",
                        "Sample_20",
                        ]
    reference = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.6m/LAN210_v0.6m.fasta"
    reference_index = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.6m/LAN210_v0.6m.idx"
    known_sites_vcf = "/run/media/mahajrod/Data/data/LAN210/fastq/210/check/reference_filtered.recode.vcf"
    run_name = "GATK"
    run_dir = "/run/media/mahajrod/Data/data/LAN210/all"
    """
    for sample_name in samples_list:
        workdir = "/run/media/mahajrod/Data/data/LAN210/all/%s/trimmed/new_alignment" % sample_name
        os.chdir(workdir)
        print("\nHandling %s...\n" % sample_name)

        add_header2bam("%s_trimmed_sorted_rm_pcr_chrom.bam" % sample_name,
                       "%s_trimmed_sorted_rm_pcr_chrom_with_header.bam" % sample_name,
                       sample_name,
                       sample_name,
                       "Illumina",
                       sample_name,
                       sample_name,
                       PICARD_dir="/home/mahajrod/Repositories/genetic/NGS_tools/picard-tools-1.115/picard-tools-1.115")

        snp_call_GATK("%s_trimmed_sorted_rm_pcr_chrom_with_header.bam" % sample_name,
                      sample_name,
                      reference,
                      known_sites_vcf,
                      stand_call_conf=100,
                      GATK_dir="/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.1-1")
        #vcftools badly parses metadata of GATK output
        os.system("vcftools --vcf %s_GATK_filtered_snps.vcf --remove-filtered-all --out %s_GATK_best_snps --recode --recode-INFO-all"
                  % (sample_name, sample_name ))
        vcf_file = "/run/media/mahajrod/Data/data/LAN210/all/%s/trimmed/new_alignment/%s_GATK_best_snps.recode.vcf" \
                   % (sample_name, sample_name)

        mutations_vcf = MutationsVcf(vcf_file, from_file=True)
        filtered_mutations, filtered_out_mutations = mutations_vcf.filter_by_reference_and_alt([("G", ["A"]), ("C", ["T"])])
        filtered_mutations.write(vcf_file[:-4] + "_filtered.vcf")
        filtered_out_mutations.write(vcf_file[:-4] + "_filtered_out.vcf")
    """
    os.chdir(run_dir)
    os.system("mkdir -p %s_vcf" % run_name)
    os.chdir("%s_vcf" % run_name)
    sub_folder_list = ["raw", "filtered_snps", "best_snps", "best_snp_only_desaminase", "best_snp_nondesaminase"]
    suffix_list = ["raw.vcf", "filtered_snps.vcf", "best_snps.recode.vcf", "best_snps.recode_filtered.vcf",	"best_snps.recode_filtered_out.vcf"]
    os.system("mkdir -p %s" % (" ".join(sub_folder_list)))
    os.chdir(run_dir)
    for subfolder, suffix in zip(sub_folder_list, suffix_list):
        os.system("cp */trimmed/new_alignment/*_GATK_%s %s_vcf/%s" % (suffix, run_name, subfolder))
