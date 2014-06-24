#/usr/bin/env python2
import os


def get_alignment(bowtie2_index,
                  sample_name,
                  min_length,
                  forward_reads,
                  forward_trim,
                  reverse_reads=None,
                  reverse_trim=None,
                  skip_correction=False):

    print("Handling %s sample..." % sample_name)
    os.system("mkdir -p trimmed")

    if reverse_reads:
        os.system("trim_galore  --length %i --phred33 --dont_gzip --clip_R1 %i --clip_R2 %i --paired -q 20 -t -o trimmed %s %s"
                  % (min_length, forward_trim, reverse_trim, forward_reads, reverse_reads))
        os.chdir("trimmed")
        os.system("fastqc -t 2 --nogroup %s_1_val_1.fq %s_2_val_2.fq" % (sample_name, sample_name))

        left_r = "spades/corrected/%s_1_val_1.00.0_0.cor.fastq" % sample_name
        right_r = "spades/corrected/%s_2_val_2.00.0_0.cor.fastq" % sample_name

        if skip_correction:
            left_r = "%s_1_val_1.fq" % sample_name
            right_r = "%s_2_val_2.fq" % sample_name
        else:
            os.system("spades.py -t 5 --only-error-correction --disable-gzip-output -1 %s_1_val_1.fq -2 %s_2_val_2.fq -o spades"
                      % (sample_name, sample_name))

        os.system("bowtie2 --phred33 -p 5 -x %s -1 %s -2 %s > %s_trimmed.sam"
                  % (bowtie2_index, left_r, right_r, sample_name))
    else:
        os.system("trim_galore  --length %i --phred33 --dont_gzip --clip_R1 %i -q 20 -t -o trimmed %s"
                  % (min_length, forward_trim, forward_reads))
        os.chdir("trimmed")
        os.system("fastqc -t 2 --nogroup %s_trimmed.fq" % (sample_name))

        unpaired_r = "spades/corrected/%s_trimmed.00.0_0.cor.fastq" % sample_name
        if skip_correction:
            unpaired_r = "%s_trimmed.fq" % sample_name
        else:
            os.system("spades.py -t 5 --only-error-correction --disable-gzip-output -s %s_trimmed.fq -o spades"
                      % sample_name)

        os.system("bowtie2 --phred33 -p 5 -x %s -U %s > %s_trimmed.sam"
                  % (bowtie2_index, unpaired_r, sample_name))

    os.system("samtools view -Sb %s_trimmed.sam | samtools sort - %s_trimmed_sorted" % (sample_name, sample_name))
    os.system("samtools rmdup %s_trimmed_sorted.bam %s_trimmed_sorted_rm_pcr.bam" % (sample_name, sample_name))
    os.system("rm -rf %s_trimmed.sam %s_trimmed_sorted.bam" % (sample_name, sample_name))
    os.system("qualimap bamqc -bam %s_trimmed_sorted_rm_pcr.bam " % sample_name)


def get_coverage_thresholds(coverage_dist_file, one_side_base_threshold=0.025, minimum_threshold=10):
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
             snp_quality=40):

    os.system("samtools mpileup  -q %i -ugf %s %s | bcftools  view -cvgN  - | vcfutils.pl varFilter -D %i -d %i > %s.vcf"
              % (alignment_quality, reference_file, alignment, max_coverage, min_coverage, sample_name))
    os.system("vcftools --vcf %s.vcf --out %s_filtered --remove-indels --recode --recode-INFO-all --minQ %i"
              % (sample_name, sample_name, snp_quality))


def snp_call_pipeline(bowtie2_index,
                      sample_name,
                      min_length,
                      reference_file,
                      right_reads,
                      right_trim,
                      left_reads=None,
                      left_trim=None,
                      skip_correction=False,
                      coverage_one_side_base_threshold=0.025,
                      coverage_minimum_threshold=10,
                      alignment_quality=40,
                      snp_quality=40):
    get_alignment(bowtie2_index, sample_name, min_length, right_reads,
                  right_trim, left_reads=left_reads, left_trim=left_trim, skip_correction=skip_correction)
    min_coverage, max_coverage = \
        get_coverage_thresholds("%s_trimmed_sorted_rm_pcr_stats/raw_data/coverage_histogram.txt" % sample_name,
                                one_side_base_threshold=coverage_one_side_base_threshold,
                                minimum_threshold=coverage_minimum_threshold)
    snp_call("%s_trimmed_sorted_rm_pcr.bam" % sample_name,
             sample_name,
             reference_file,
             min_coverage,
             max_coverage,
             alignment_quality=alignment_quality,
             snp_quality=snp_quality)