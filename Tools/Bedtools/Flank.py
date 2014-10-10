#!/usr/bin/env python
from Tools.Abstract import Tool


class Flank(Tool):
    """
    Tool:    bedtools flank (aka flankBed)
    Version: v2.17.0
    Summary: Creates flanking interval(s) for each BED/GFF/VCF feature.

    Usage:   bedtools flank [OPTIONS] -i <bed/gff/vcf> -g <genome> [-b <int> or (-l and -r)]

    Options:
        -b	Create flanking interval(s) using -b base pairs in each direction.
            - (Integer) or (Float, e.g. 0.1) if used with -pct.

        -l	The number of base pairs that a flank should start from
            orig. start coordinate.
            - (Integer) or (Float, e.g. 0.1) if used with -pct.

        -r	The number of base pairs that a flank should end from
            orig. end coordinate.
            - (Integer) or (Float, e.g. 0.1) if used with -pct.

        -s	Define -l and -r based on strand.
            E.g. if used, -l 500 for a negative-stranded feature,
            it will start the flank 500 bp downstream.  Default = false.

        -pct	Define -l and -r as a fraction of the feature's length.
            E.g. if used on a 1000bp feature, -l 0.50,
            will add 500 bp "upstream".  Default = false.

        -header	Print the header from the input file prior to results.

    Notes:
        (1)  Starts will be set to 0 if options would force it below 0.
        (2)  Ends will be set to the chromosome length if requested flank would
        force it above the max chrom length.
        (3)  In contrast to slop, which _extends_ intervals, bedtools flank
        creates new intervals from the regions just up- and down-stream
        of your existing intervals.
        (4)  The genome file should tab delimited and structured as follows:

        <chromName><TAB><chromSize>

        For example, Human (hg19):
        chr1	249250621
        chr2	243199373
        ...
        chr18_gl000207_random	4262

    Tips:
        One can use the UCSC Genome Browser's MySQL database to extract
        chromosome sizes. For example, H. sapiens:

        mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -e \
        "select chrom, size from hg19.chromInfo"  > hg19.genome
    """

    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "bedtools flank", path=path, max_threads=max_threads)

    def get(self, in_file, genome_file, left, right, fraction_mode=False, strand_based=False,
            print_header=False, out_file=None):

        options = ""
        options += " -header" if print_header else ""
        options += " -s" if strand_based else ""
        options += " -pct" if fraction_mode else ""
        options += " -l %s" % str(left)
        options += " -r %s" % str(right)
        options += " -i %s" % in_file
        options += " -g %s" % genome_file
        options += " > %s" % out_file if out_file else ""

        self.execute(options)

