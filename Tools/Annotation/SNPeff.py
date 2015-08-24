#!/usr/bin/env python

from Tools.Abstract import JavaTool
from collections import Iterable
from Routines.Functions import check_path


class SNPeff(JavaTool):
    def __init__(self, max_threads=4, jar_path="", jar="snpEff.jar",
                 max_memory="500m", config_path=""):
        JavaTool.__init__(self, jar, java_path="", max_threads=max_threads, jar_path=jar_path, max_memory=max_memory)
        self.config_path = config_path


    def build_db(self, datafiles, input_format="gff3", build_type=None):
        """
        To build custom db at first add record to config file

        # Mouse genome, version mm37.61
        mm37.61.genome : Mouse

        Build DB options:
            input_format:
            -embl                   : Use Embl format. It implies '-1'.
            -genbank                : Use GenBank format. It implies '-1'.
            -gff2                   : Use GFF2 format (obsolete). It implies '-1'.
            -gff3                   : Use GFF3 format. It implies '-1'
            -gtf22                  : Use GTF 2.2 format. It implies '-1'. Default
            -refseq                 : Use RefSeq table from UCSC. It implies '-0'.
            -txt                    : Use TXT format (obsolete).

            build_type
            -knowngenes             : Use KnownGenes table from UCSC. It implies '-0'.
            -onlyReg                : Only build regulation tracks.
            -cellType <type>        : Only build regulation tracks for cellType <type>.

        Generic options:
            -0                      : File positions are zero-based (same as '-inOffset 0 -outOffset 0')
            -1                      : File positions are one-based (same as '-inOffset 1 -outOffset 1')
            -c , -config            : Specify config file
            -h , -help              : Show this help and exit
            -if, -inOffset          : Offset input by a number of bases. E.g. '-inOffset 1' for one-based input files
            -of, -outOffset         : Offset output by a number of bases. E.g. '-outOffset 1' for one-based output files
            -noLog                  : Do not report usage statistics to server
            -q , -quiet             : Quiet mode (do not show any messages or errors)
            -v , -verbose           : Verbose mode
        """

        options = "-v -%s" % input_format
        options += " %s" % build_type if build_type else ""
        options += " " + " ".join(datafiles) if isinstance(datafiles, list) else " " + datafiles

        self.execute(options, cmd="build")

    def annotate(self, genome, input_vcf, output_vcf, summary_file=None, verbose=True, no_downstream=False,
                 no_intergenic=False, no_intron=False, no_upstream=False, no_utr=False):

        options = " -v" if verbose else ""
        options += " -s %s" % summary_file if summary_file else ""
        options += " -no-downstream" if no_downstream else ""
        options += " -no-intergenic" if no_intergenic else ""
        options += " -no-intron" if no_intron else ""
        options += " -no-upstream" if no_upstream else ""
        options += " -no-utr" if no_utr else ""
        options += " %s" % genome
        options += " %s" % input_vcf
        options += " > %s" % output_vcf

        self.execute(options, cmd="")


if __name__ == "__main__":
    import os
    workdir = "/home/mahajrod/genetics/desaminases/data/S288C_R64/SNPeff"
    SNPeff_path = "/home/mahajrod/Repositories/genetic/NGS_tools/snpEff"
    SNPeff = SNPeff(jar_path=SNPeff_path)

    #print(SNPeff.jar_path)
    #os.chdir(workdir)
    #gff_file = "/home/mahajrod/genetics/desaminases/data/S288C_R64/merged_saccharomyces_cerevisiae_R64-1-1_20110208_Nagalakshmi_UTR_2008.gff3"
    #SNPeff.build_db(gff_file, input_format="gff3")

    """
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/intersect_S288C/normal_regions"
    ref_genome = "S288C_R64_UTR"
    in_file = "in_normal_regions_homo.vcf"
    out_file_suffix = "_S288C_R64_UTR_annotated.vcf"
    """

    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/fastq/intersect_LAN210_v0.10m/normal_regions"
    ref_genome = "LAN210_v0.10m"
    in_file = "variants_normal_regions_homo.vcf"
    out_file_suffix = "_LAN210_v0.10m_annotated.vcf"

    os.chdir(workdir)

    #os.system("mkdir -p %s" % ref_genome)
    #os.chdir(ref_genome)

    SNPeff.annotate(ref_genome, "../" + in_file, in_file[:-4] + out_file_suffix)