#!/usr/bin/env python

import os
from Tools.Abstract import Tool
from Routines.Functions import check_path


class Intersect(Tool):
    """
    Tool:    bedtools intersect (aka intersectBed)
    Version: v2.17.0
    Summary: Report overlaps between two feature files.

    Usage:   bedtools intersect [OPTIONS] -a <bed/gff/vcf> -b <bed/gff/vcf>

    Options:
        -abam	The A input file is in BAM format.  Output will be BAM as well.

        -ubam	Write uncompressed BAM output. Default writes compressed BAM.

        -bed	When using BAM input (-abam), write output as BED. The default
            is to write output in BAM when using -abam.

        -wa	Write the original entry in A for each overlap.

        -wb	Write the original entry in B for each overlap.
            - Useful for knowing _what_ A overlaps. Restricted by -f and -r.

        -loj	Perform a "left outer join". That is, for each feature in A
            report each overlap with B.  If no overlaps are found,
            report a NULL feature for B.

        -wo	Write the original A and B entries plus the number of base
            pairs of overlap between the two features.
            - Overlaps restricted by -f and -r.
              Only A features with overlap are reported.

        -wao	Write the original A and B entries plus the number of base
            pairs of overlap between the two features.
            - Overlapping features restricted by -f and -r.
              However, A features w/o overlap are also reported
              with a NULL B feature and overlap = 0.

        -u	Write the original A entry _once_ if _any_ overlaps found in B.
            - In other words, just report the fact >=1 hit was found.
            - Overlaps restricted by -f and -r.

        -c	For each entry in A, report the number of overlaps with B.
            - Reports 0 for A entries that have no overlap with B.
            - Overlaps restricted by -f and -r.

        -v	Only report those entries in A that have _no overlaps_ with B.
            - Similar to "grep -v" (an homage).

        -f	Minimum overlap required as a fraction of A.
            - Default is 1E-9 (i.e., 1bp).
            - FLOAT (e.g. 0.50)

        -r	Require that the fraction overlap be reciprocal for A and B.
            - In other words, if -f is 0.90 and -r is used, this requires
              that B overlap 90% of A and A _also_ overlaps 90% of B.

        -s	Require same strandedness.  That is, only report hits in B
            that overlap A on the _same_ strand.
            - By default, overlaps are reported without respect to strand.

        -S	Require different strandedness.  That is, only report hits in B
            that overlap A on the _opposite_ strand.
            - By default, overlaps are reported without respect to strand.

        -split	Treat "split" BAM or BED12 entries as distinct BED intervals.

        -sorted	Use the "chromsweep" algorithm for sorted (-k1,1 -k2,2n) input

        -header	Print the header from the A file prior to results.

    Notes:
        (1) When a BAM file is used for the A file, the alignment is retained if overlaps exist,
        and exlcuded if an overlap cannot be found.  If multiple overlaps exist, they are not
        reported, as we are only testing for one or more overlaps.
    """

    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "bedtools intersect", path=path, max_threads=max_threads)

    def intersect(self, a_file, b_file, out_file, method="-v"):

        options = " -a %s" % a_file
        options += " -b %s" % b_file
        options += " %s" % method
        options += " > %s" % out_file

        self.execute(options)

    def intersect_set(self, file_list, outfile, method="-v"):
        os.system("cat %s > temp0" % (file_list[0]))
        for i in range(1, len(file_list)):
            os.system("%sbedtools intersect %s -a %s -b %s > %s"
                      % (self.bedtools_dir, method, "temp%i" % (i-1), file_list[i], "temp%i" % i))
        os.system("cat temp%i > %s" % (i, outfile))
        os.system("rm -f %s" % " ".join(["temp%i" % x for x in range(0, i+1)]))
