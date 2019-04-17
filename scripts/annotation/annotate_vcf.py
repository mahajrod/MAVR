#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Annotation import SNPeff

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_vcf", action="store", dest="input", required=True,
                    help="Input .vcf file")
parser.add_argument("-o", "--output_vcf", action="store", dest="output", required=True,
                    help="Annotated output .vcf file")
parser.add_argument("-g", "--genome", action="store", dest="genome", required=True,
                    help="Genome to use (from SNPeff database)")
parser.add_argument("-d", "--snpeff_dir", action="store", dest="snpeff_dir",
                    default="/home/mahajrod/Repositories/genetic/NGS_tools/snpEff/",
                    help="Directory with SNPeff jar")
parser.add_argument("-s", "--summary_file", action="store", dest="summary_file",
                    help="Name of stats file (summary). Default is 'snpEff_summary.html'")
parser.add_argument("--no_downstream", action="store_true", dest="no_downstream",
                    help="Do not show DOWNSTREAM changes")
parser.add_argument("--no_intergenic", action="store_true", dest="no_intergenic",
                    help="Do not show INTERGENIC changes")
parser.add_argument("--no_intron", action="store_true", dest="no_intron",
                    help="Do not show INTRON changes")
parser.add_argument("--no_upstream", action="store_true", dest="no_upstream",
                    help="Do not show UPSTREAM changes")
parser.add_argument("--no_utr", action="store_true", dest="no_utr",
                    help="Do not show 5_PRIME_UTR or 3_PRIME_UTR changes")
parser.add_argument("-v", "--verbose", action="store_true", dest="verbose",
                    help="Verbose mode")
parser.add_argument("-m", "--memory", action="store", dest="memory",
                    help="Memory limit for java. Default: 500m")

args = parser.parse_args()

SNPeff.jar_path = args.snpeff_dir
SNPeff.max_memory = args.memory
SNPeff.annotate(args.genome, args.input, args.output, summary_file=args.summary_file, verbose=args.verbose,
                no_downstream=args.no_downstream, no_intergenic=args.no_intergenic, no_intron=args.no_intron,
                no_upstream=args.no_upstream, no_utr=args.no_utr)
