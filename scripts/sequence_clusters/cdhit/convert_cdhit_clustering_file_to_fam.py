#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Tools.Clustering import CDHit


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--clustering_file", action="store", dest="clustering_file", required=True,
                    help="Input file with CDHit clusters")
parser.add_argument("-o", "--output", action="store", dest="output", required=True,
                    help="Output fam file")


args = parser.parse_args()


CDHit.convert_clustering_to_fam(args.clustering_file, args.output)
