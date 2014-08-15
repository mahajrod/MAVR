#/usr/bin/env python
import os
import collections
from Bio import SeqIO
from General import check_path
from Tools.Picard import add_header2bam
#from MutAnalysis.Mutation import *
from Tools.AssemblyTools import spades
from Tools.FilterTools import trim_galore
from Tools.AlignmentTools import Bowtie2, BWA


# TODO: refactor, remove absolete functions and so on


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
                                    remove_nonPA=False):
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
    os.system("qualimap bamqc -bam %s_trimmed_sorted_rm_pcr_chrom.bam " % sample_name)

    if mitochondrial_bed_file:
        os.system("samtools view -b -F %i -L %s %s_trimmed_sorted_rm_pcr.bam  > %s_trimmed_sorted_rm_pcr_mt.bam"
              % (filter_parameter, mitochondrial_bed_file, sample_name, sample_name))
        os.system("samtools index %s_trimmed_sorted_rm_pcr_mt.bam" % sample_name)
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
                  find_separated_alignments=True):

    print("Handling %s sample..." % sample_name)
    output_file = "%s_trimmed.sam" % sample_name
    if aligner == "bowtie2":
        bowtie2 = Bowtie2()
        bowtie2.align(index, forward_reads, reverse_reads,
                      max_threads=max_threads, quality_score=quality_score,
                      alignment_mode="very-sensitive", output_file=output_file,
                      find_discordant_alignments=find_discordant_alignments,
                      find_separated_alignments=find_separated_alignments)
    elif aligner == "bwa-mem":
        bwa = BWA()
        bwa.align_mem(index, forward_reads, reverse_reads, output_file=output_file, max_threads=max_threads)
    """
    if reverse_reads:
        os.system("bowtie2 --very-sensitive --%s -p %i -x %s -1 %s -2 %s > %s_trimmed.sam"
                  % (quality_score, max_threads, bowtie2_index, forward_reads, reverse_reads, sample_name))
    else:
        os.system("bowtie2 --very-sensitive --%s -p %i -x %s -U %s > %s_trimmed.sam"
                  % (quality_score, max_threads, bowtie2_index, forward_reads, sample_name))
    """
    alignment_sorting_and_filtering(sample_name, output_file, chromosomes_bed_file, mitochondrial_bed_file,
                                    remove_SA=remove_SA, remove_nonPA=remove_nonPA)


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
"""
def correct_reads(forward_reads,
                  reverse_reads=None,
                  max_threads=4,
                  ion_torrent=False):
    #print("Handling %s sample..." % sample_name)
    os.system("mkdir -p spades")
    platform = ""
    if ion_torrent:
        platform = "--iontorrent"
    if reverse_reads:
        os.system("spades.py -t %i %s --only-error-correction --disable-gzip-output -1 %s -2 %s -o spades"
                  % (max_threads, platform, forward_reads, reverse_reads))
    else:
        os.system("spades.py -t %i %s --only-error-correction --disable-gzip-output -s %s -o spades"
                  % (max_threads, platform, forward_reads))
    os.chdir("spades/corrected")
    os.system("fastqc -t 2 --nogroup *.f*q")

def assembly_reads(forward_reads,
                  reverse_reads=None,
                  max_threads=4,
                  ion_torrent=False):
    #print("Handling %s sample..." % sample_name)
    os.system("mkdir -p spades")
    platform = ""
    if ion_torrent:
        platform = "--iontorrent"
    if reverse_reads:
        os.system("spades.py -t %i %s --disable-gzip-output -1 %s -2 %s -o spades"
                  % (max_threads, platform, forward_reads, reverse_reads))
    else:
        os.system("spades.py -t %i %s --disable-gzip-output -s %s -o spades"
                  % (max_threads, platform, forward_reads))
    os.chdir("spades/corrected")
    os.system("fastqc -t 2 --nogroup *.f*q")

def filter_data(sample_name,
                min_length,
                forward_reads,
                forward_trim,
                reverse_reads=None,
                reverse_trim=None,
                quality_score="phred33",
                adapter="AGATCGGAAGAGC",
                quality_treshold=20):
    print("Handling %s sample..." % sample_name)
    os.system("mkdir -p trimmed")
    # if forward_trim == None skip 5' trimming
    trim_str = ""

    if reverse_reads:
        if forward_trim:
            trim_str = "--clip_R1 %i --clip_R2 %i" % (forward_trim, reverse_trim)
        os.system("trim_galore -a %s  --length %i --%s --dont_gzip %s --paired -q %i -o trimmed %s %s"
                  % (adapter, min_length, quality_score, trim_str, quality_treshold, forward_reads, reverse_reads))
    else:
        if forward_trim:
            trim_str = "--clip_R1 %i" % forward_trim
        os.system("trim_galore -a %s  --length %i --%s --dont_gzip %s -q %i -o trimmed %s"
                  % (adapter, min_length, quality_score, trim_str, quality_treshold, forward_reads))
    os.chdir("trimmed")
    os.system("fastqc -t 2 --nogroup *.fq")
"""


def correct_reads(forward_reads,
                  reverse_reads=None,
                  max_threads=4,
                  platform="illumina",
                  output_dir="spades"):
    #print("Handling %s sample..." % sample_name)
    spades(forward_reads, reverse_reads=reverse_reads, only_correction=True,
           max_threads=max_threads, platform=platform,
           gzip_corrected_reads=False, output_dir=output_dir)


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

"""
def get_alignment(bowtie2_index,
                  sample_name,
                  min_length,
                  forward_reads,
                  forward_trim,
                  chromosomes_bed_file,
                  mitochondrial_bed_file,
                  reverse_reads=None,
                  reverse_trim=None,
                  skip_correction=False,
                  max_threads=5,
                  quality_score="phred33",
                  adapter="AGATCGGAAGAGC"):
    #Illumina standard adapter AGATCGGAAGAGC
    #Nextera adapter CTGTCTCTTATACACATCT
    #Ion torrent adapter P1 CCACTACGCCTCCGC
    print("Handling %s sample..." % sample_name)
    os.system("mkdir -p trimmed")

    if reverse_reads:
        os.system("trim_galore -a %s  --length %i --%s --dont_gzip --clip_R1 %i --clip_R2 %i --paired -q 20 -t -o trimmed %s %s"
                  % (adapter, min_length, quality_score, forward_trim, reverse_trim, forward_reads, reverse_reads))
        os.chdir("trimmed")
        os.system("fastqc -t 2 --nogroup %s_1_val_1.fq %s_2_val_2.fq" % (sample_name, sample_name))

        left_r = "spades/corrected/%s_1_val_1.00.0_0.cor.fastq" % sample_name
        right_r = "spades/corrected/%s_2_val_2.00.0_0.cor.fastq" % sample_name

        if skip_correction:
            left_r = "%s_1_val_1.fq" % sample_name
            right_r = "%s_2_val_2.fq" % sample_name
        else:
            os.system("spades.py -t %i --only-error-correction --disable-gzip-output -1 %s_1_val_1.fq -2 %s_2_val_2.fq -o spades"
                      % (max_threads,sample_name, sample_name))

        os.system("bowtie2 --%s -p %i -x %s -1 %s -2 %s > %s_trimmed.sam"
                  % (quality_score, max_threads, bowtie2_index, left_r, right_r, sample_name))
    else:
        os.system("trim_galore  --length %i --%s --dont_gzip --clip_R1 %i -q 20 -t -o trimmed %s"
                  % (min_length, quality_score, forward_trim, forward_reads))
        os.chdir("trimmed")
        os.system("fastqc -t 2 --nogroup %s_trimmed.fq" % (sample_name))

        unpaired_r = "spades/corrected/%s_trimmed.00.0_0.cor.fastq" % sample_name
        if skip_correction:
            unpaired_r = "%s_trimmed.fq" % sample_name
        else:
            os.system("spades.py -t %i --only-error-correction --disable-gzip-output -s %s_trimmed.fq -o spades"
                      % (max_threads, sample_name))

        os.system("bowtie2 --%s -p %i -x %s -U %s > %s_trimmed.sam"
                  % (quality_score, max_threads, bowtie2_index, unpaired_r, sample_name))

    alignment_sorting_and_filtering(sample_name, chromosomes_bed_file, mitochondrial_bed_file)
"""

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
    gatk_dir = check_path(GATK_dir)
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
