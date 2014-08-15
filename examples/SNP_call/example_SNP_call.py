#!/usr/bin/env python
from Pipelines.SNPCall import *
from Tools.Picard import *
from Tools.GATK import UnifiedGenotyper, SelectVariants, VariantFiltration, CombineVariants
from Parser.VCF import CollectionVCF
from Routines.Functions import check_path
#in case reads are already trimmed and corrected


def get_read_filenames(sample_name, read_type, sample_dir="."):
    dir_list = []
    for file_name in os.listdir(sample_dir):
        if file_name[-3:] == ".fq" or file_name[-6:] == ".fastq":
            dir_list.append(file_name)
    if read_type == "SE":
        for filename in dir_list:
            if ("_R1" in filename or ".R1" in filename) and (sample_name in filename):
                return filename, None
    for filename in dir_list:
        if ("_R1" in filename) and (sample_name in filename):
            left_reads = filename
        if ("_R2" in filename) and (sample_name in filename):
            right_reads = filename
    return left_reads, right_reads

"""
reference_dir = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.8m/"
reference = "%sLAN210_v0.8m.fasta" % reference_dir
reference_dict = "%sLAN210_v0.8m.dict" % reference_dir
reference_index = "%sLAN210_v0.8m" % reference_dir
"""
"""
reference_name = "LAN210_v0.9m"
reference_dir = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/"
"""
reference_name = "LAN210_v0.10m"
reference_dir = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/"
reference = "%s%s.fasta" % (reference_dir, reference_name)
reference_dict = "%s%s.dict" % (reference_dir, reference_name)
reference_index = "%s%s" % (reference_dir, reference_name)

workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all"

gatk_dir = check_path("/home/mahajrod/Repositories/genetic/NGS_tools/GenomeAnalysisTK-3.2-0/")
platform = "illumina"
read_subdir = "trimmed/spades/corrected/"

#samples_file = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/samples_nonwt.t"
#samples_file = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/samples_wt.t"
samples_file = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/samples_nonHAP.t"
#alignment_dir = "alignment"
#alignment_dir = "alignment_LAN210_v0.9m"
alignment_dir = "alignment_LAN210_v0.10m"
os.chdir(reference_dir)
make_fasta_dict(reference, reference_dict, PICARD_dir="/home/mahajrod/Repositories/genetic/NGS_tools/picard-tools-1.115/picard-tools-1.115/")
get_chromosomes_bed(reference, reference_index, mitochondrial_region_name="mt",
                        chrom_out_file="chromosomes.bed", mito_out_file="mt.bed", reference_filetype="fasta")
os.system("samtools faidx %s" % reference)

os.chdir(workdir)


samples_list = []
with open(samples_file, "r") as in_fd:
    in_fd.readline()
    for line in in_fd:
        if line == "\n":
            continue
        tmp = line.strip().split("\t")
        tmp[4] = int(tmp[4])
        tmp[5] = int(tmp[5])
        if tmp[6] == ".":
            tmp[6] = None
        else:
            tmp[6] = int(tmp[6])
        if tmp[7] == ".":
            tmp[7] = None
        else:
            tmp[7] = int(tmp[7])
        #N001-LAN200-wt-NA-NA-RUN1    illumina    TruSeq    PE    101    5234888    6    6    AGATCGGAAGAGC
        samples_list.append(tmp)

for sample_name, platform, library, read_type, read_length, n_of_reads, left_trim, right_trim, adapter in samples_list:
    os.chdir(workdir)
    os.chdir(sample_name)
    os.system("mkdir -p %s" % alignment_dir)
    left_reads, right_reads = get_read_filenames(sample_name, read_type, sample_dir=read_subdir)
    forward_reads = "../" + read_subdir + left_reads
    reverse_reads = None
    if right_reads:
        reverse_reads = "../" + read_subdir + right_reads
    os.chdir(alignment_dir)
    get_alignment(reference_index,
                  sample_name,
                  forward_reads,
                  "%schromosomes.bed" % reference_dir,
                  "%smt.bed" % reference_dir,
                  reverse_reads=reverse_reads,
                  max_threads=5,
                  quality_score="phred33",
                  remove_SA=True,
                  remove_nonPA=True)

    add_header2bam("%s_trimmed_sorted_rm_pcr_chrom.bam" % sample_name,
                "%s_trimmed_sorted_rm_pcr_chrom_with_header.bam" % sample_name,
                sample_name,
                sample_name,
                platform,
                sample_name,
                sample_name,
                PICARD_dir="/home/mahajrod/Repositories/genetic/NGS_tools/picard-tools-1.115/picard-tools-1.115")

    vcf_all = "%s_GATK_raw.vcf" % sample_name
    vcf_indel = "%s_GATK_raw_indel.vcf" % sample_name
    vcf_SNP = "%s_GATK_raw_SNP.vcf" % sample_name
    vcf_filtered_indel = "%s_GATK_filtered_indel.vcf" % sample_name
    vcf_filtered_SNP = "%s_GATK_filtered_SNP.vcf" % sample_name
    vcf_best_indel = "%s_GATK_best_indel.vcf" % sample_name
    vcf_best_SNP = "%s_GATK_best_SNP.vcf" % sample_name
    vcf_best_merged = "%s_GATK_best_merged.vcf" % sample_name
    vcf_best_merged_hetero = "%s_GATK_best_merged_hetero.vcf" % sample_name
    vcf_best_merged_homo = "%s_GATK_best_merged_homo.vcf" % sample_name

    UnifiedGenotyper.variant_call("%s_trimmed_sorted_rm_pcr_chrom_with_header.bam" % sample_name,
                                  reference,
                                  stand_emit_conf=40,
                                  stand_call_conf=100,
                                  GATK_dir=gatk_dir,
                                  num_of_threads=5,
                                  output_mode="EMIT_VARIANTS_ONLY",
                                  discovery_mode="BOTH",
                                  output_file=vcf_all)

    SelectVariants.get_indel(gatk_dir, reference, vcf_all, vcf_indel)
    SelectVariants.get_SNP(gatk_dir, reference, vcf_all, vcf_SNP)

    VariantFiltration.filter_bad_SNP(gatk_dir, reference, vcf_SNP, vcf_filtered_SNP)
    VariantFiltration.filter_bad_indel(gatk_dir, reference, vcf_indel, vcf_filtered_indel)
    SelectVariants.remove_filtered(gatk_dir, reference, vcf_filtered_SNP, vcf_best_SNP)
    SelectVariants.remove_filtered(gatk_dir, reference, vcf_filtered_indel, vcf_best_indel)

    CombineVariants.combine_from_same_source(gatk_dir, reference, [vcf_best_SNP, vcf_best_indel], vcf_best_merged)

    best_merged = CollectionVCF(vcf_file=vcf_best_merged)
    best_merged_homo, best_merged_hetero = best_merged.split_by_zygoty()
    best_merged_homo.write(vcf_best_merged_homo)
    best_merged_hetero.write(vcf_best_merged_hetero)




