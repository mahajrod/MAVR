#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import argparse
from RouToolPa.Routines import NCBIRoutines

parser = argparse.ArgumentParser()

parser.add_argument("-t", "--taxa", action="store", dest="taxa", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of taxon names for query")
parser.add_argument("-e", "--email", action="store", dest="email", required=True,
                    help="Email used in Entrez queues")
parser.add_argument("-d", "--output_dir", action="store", dest="output_dir", default="./",
                    help="Directory to write output. Default: current directory")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", default="",
                    help="Prefix of output files. Default: no prefix")
parser.add_argument("--min_scaffold_n50", action="store", dest="min_scaffold_n50", type=int,
                    help="Minimum N50 of assembly scaffolds to pass filters. Default: not set")
parser.add_argument("--min_contig_n50", action="store", dest="min_contig_n50", type=int,
                    help="Minimum N50 of assembly contiigs to pass filters. Default: not set")
parser.add_argument("--max_scaffold_l50", action="store", dest="max_scaffold_l50", type=int,
                    help="Maximum L50 of assembly scaffolds to pass filters. Default: not set")
parser.add_argument("--max_contig_l50", action="store", dest="max_contig_l50", type=int,
                    help="Maximum L50 of assembly contigs to pass filters. Default: not set")
parser.add_argument("--max_contig_count", action="store", dest="max_contig_count", type=int,
                    help="Maximum number of assembly contigs to pass filters. Default: not set")
parser.add_argument("--max_scaffold_count", action="store", dest="max_scaffold_count", type=int,
                    help="Maximum number of assembly scaffolds to pass filters. Default: not set")
parser.add_argument("--max_chromosome_count", action="store", dest="max_chromosome_count", type=int,
                    help="Maximum number of assembly chromosomes to pass filters. Default: not set")
parser.add_argument("--min_chromosome_count", action="store", dest="min_chromosome_count", type=int,
                    help="Minimum number of assembly chromosomes to pass filters. Default: not set")
parser.add_argument("--max_unlocalized_scaffolds", action="store", dest="max_unlocalized_scaffolds", type=int,
                    help="Maximum number of unlocalized scaffolds in assembly to pass filters. Default: not set")
parser.add_argument("--max_unplaced_scaffolds", action="store", dest="max_unplaced_scaffolds", type=int,
                    help="Maximum number of unplaced scaffolds in assembly to pass filters. Default: not set")
parser.add_argument("--max_total_length", action="store", dest="max_total_length", type=int,
                    help="Maximum total length of assembly to pass filters. Default: not set")
parser.add_argument("--min_total_length", action="store", dest="min_total_length", type=int,
                    help="Minimum total length of assembly to pass filters. Default: not set")
parser.add_argument("--max_ungapped_length", action="store", dest="max_ungapped_length", type=int,
                    help="Maximum ungapped length of assembly to pass filters. Default: not set")
parser.add_argument("--min_ungapped_length", action="store", dest="min_ungapped_length", type=int,
                    help="Minimum ungapped length of assembly to pass filters. Default: not set")
parser.add_argument("--include_ambiguous_species", action="store_true", dest="include_ambiguous_species", default=False,
                    help="Include ambiguous species (hybrids, candidates, unknown and not identified to species level)."
                         " Default: false")

args = parser.parse_args()

NCBIRoutines.get_taxa_genomes_summary(args.taxa,
                                      args.email,
                                      args.output_dir,
                                      args.output_prefix,
                                      min_scaffold_n50=args.min_scaffold_n50,
                                      min_contig_n50=args.min_contig_n50,
                                      max_scaffold_l50=args.max_scaffold_l50,
                                      max_contig_l50=args.max_contig_l50,
                                      max_contig_count=args.max_contig_count,
                                      max_scaffold_count=args.max_scaffold_count,
                                      max_chromosome_count=args.max_chromosome_count,
                                      min_chromosome_count=args.min_chromosome_count,
                                      max_unlocalized_scaffolds=args.max_unlocalized_scaffolds,
                                      max_unplaced_scaffolds=args.max_unplaced_scaffolds,
                                      max_total_length=args.max_total_length,
                                      min_total_length=args.min_total_length,
                                      max_ungapped_length=args.max_ungapped_length,
                                      min_ungapped_length=args.min_ungapped_length,
                                      no_ambiguous_species=not args.include_ambiguous_species)
