#!/usr/bin/env python
import os

#TODO: rewrite as package
def spades(forward_reads,
           reverse_reads=None,
           sanger_reads=None,
           pacbio_reads=None,
           only_correction=False,
           max_threads=4,
           platform="illumina",
           gzip_corrected_reads=True,
           output_dir="spades",
           kmer_length_list=[]):
    #suitable for spades v3.0.0
    """
    SPAdes genome assembler v.3.0.0

    Usage: spades.py [options] -o <output_dir>

    Basic options:
    -o	<output_dir>	directory to store all the resulting files (required)
    --sc			this flag is required for MDA (single-cell) data
    --iontorrent		this flag is required for IonTorrent data
    --test			runs SPAdes on toy dataset
    -h/--help		prints this usage message

    Input data:
    --12	<filename>	file with interlaced forward and reverse paired-end reads
    -1	<filename>	file with forward paired-end reads
    -2	<filename>	file with reverse paired-end reads
    -s	<filename>	file with unpaired reads
    --pe<#>-12	<filename>	file with interlaced reads for paired-end library number <#> (<#> = 1,2,3,4,5)
    --pe<#>-1	<filename>	file with forward reads for paired-end library number <#> (<#> = 1,2,3,4,5)
    --pe<#>-2	<filename>	file with reverse reads for paired-end library number <#> (<#> = 1,2,3,4,5)
    --pe<#>-s	<filename>	file with unpaired reads for paired-end library number <#> (<#> = 1,2,3,4,5)
    --pe<#>-<or>	orientation of reads for paired-end library number <#> (<#> = 1,2,3,4,5; <or> = fr, rf, ff)
    --mp<#>-12	<filename>	file with interlaced reads for mate-pair library number <#> (<#> = 1,2,3,4,5)
    --mp<#>-1	<filename>	file with forward reads for mate-pair library number <#> (<#> = 1,2,3,4,5)
    --mp<#>-2	<filename>	file with reverse reads for mate-pair library number <#> (<#> = 1,2,3,4,5)
    --mp<#>-s	<filename>	file with unpaired reads for mate-pair library number <#> (<#> = 1,2,3,4,5)
    --mp<#>-<or>	orientation of reads for mate-pair library number <#> (<#> = 1,2,3,4,5; <or> = fr, rf, ff)
    --sanger	<filename>	file with Sanger reads
    --pacbio	<filename>	file with PacBio reads
    --trusted-contigs	<filename>	file with trusted contigs
    --untrusted-contigs	<filename>	file with untrusted contigs

    Pipeline options:
    --only-error-correction	runs only read error correction (without assembling)
    --only-assembler	runs only assembling (without read error correction)
    --careful		tries to reduce number of mismatches and short indels
    --continue		continue run from the last available check-point
    --restart-from	<cp>	restart run with updated options and from the specified check-point ('ec', 'as', 'k<int>', 'mc')
    --disable-gzip-output	forces error correction not to compress the corrected reads
    --disable-rr		disables repeat resolution stage of assembling

    Advanced options:
    --dataset	<filename>	file with dataset description in YAML format
    -t/--threads	<int>		number of threads
                    [default: 16]
    -m/--memory	<int>		RAM limit for SPAdes in Gb (terminates if exceeded)
                    [default: 250]
    --tmp-dir	<dirname>	directory for temporary files
                    [default: <output_dir>/tmp]
    -k		<int,int,...>	comma-separated list of k-mer sizes (must be odd and
                    less than 128) [default: 'auto']
    --phred-offset	<33 or 64>	PHRED quality offset in the input reads (33 or 64)
                    [default: auto-detect]
    """

    os.system("mkdir -p %s" % output_dir)

    if platform == "illumina":
        platform = ""
    elif platform == "ion_torrent":
        platform = "--iontorrent"

    if reverse_reads:
        reads = "-1 %s -2 %s" % (forward_reads, reverse_reads)
    else:
        reads = "-s %s" % forward_reads

    if sanger_reads:
        reads += " --sanger %s" % sanger_reads
    if pacbio_reads:
        reads += " --pacbio %s" % pacbio_reads

    if only_correction:
        mode = " --only-error-correction"
    else:
        mode = ""

    if gzip_corrected_reads:
        gzip_output = ""
    else:
        gzip_output = "--disable-gzip-output"

    kmers = ""
    if kmer_length_list:
        kmers = "-k " + ",".join(map(lambda x: str(x), kmer_length_list))

    os.system("spades.py -t %i %s %s %s %s --careful %s -o %s"
              % (max_threads, mode, platform, gzip_output, reads, kmers, output_dir))

    if not gzip_corrected_reads:
        os.chdir("%s/corrected" % output_dir)
        os.system("fastqc -t 2 --nogroup *.f*q")
        os.system("cd ../..")