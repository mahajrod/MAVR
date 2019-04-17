#!/usr/bin/env python
__author__ = 'mahajrod'
import os
import argparse
import pprint
from BCBio.GFF import GFFExaminer
from BCBio import GFF
import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
plt.ioff()

from RouToolPa.Parsers.VCF import CollectionVCF
parser = argparse.ArgumentParser()

parser.add_argument("-f", "--gff", action="store", dest="gff",
                    help="gff with annotations")
parser.add_argument("-v", "--vcf", action="store", dest="vcf",
                    help="vcf with variants")
parser.add_argument("-r", "--right_region", action="store", dest="right",
                    help="size of region right to CDS start", type=int, default=0)
parser.add_argument("-l", "--left_region", action="store", dest="left",
                    help="size of region left to CDS start", type=int, default=0)
parser.add_argument("-b", "--both_region", action="store", dest="both",
                    help="size of regions left and right to CDS start< overrides -l and -r options",
                    type=int, default=None)
parser.add_argument("-g", "--gene_variants", action="store", dest="gene_variants",
                    help="Positions of variants for each gene",
                    default="gene_variants_positions.t")

parser.add_argument("-s", "--start_histogramm_prefix", action="store", dest="start_histogramm_prefix",
                    help="Start histogramm file prefix",
                    default="start_histogramm")

parser.add_argument("-e", "--end_histogramm_prefix", action="store", dest="end_histogramm_prefix",
                    help="End histogramm file prefix",
                    default="end_histogramm")

args = parser.parse_args()

if args.both is not None:
    args.left = args.both
    args.right = args.both

if (not args.right) and (not args.left):
    raise ValueError("Both left and right regions were not set")

with open(args.gff, "r") as in_fd:
    record_dict = dict([(record.id, record) for record in GFF.parse(in_fd)])

variants = CollectionVCF(from_file=True, vcf_file=args.vcf)
gene_variants_positions = []
all_variant_start_positions = []
all_variant_end_positions = []
print(args.left)
print(args.right)
for record_id in record_dict:
    for feature in record_dict[record_id].features:
        if feature.type != "gene":
            continue
        #print(feature.sub_features)
        for sub_feature in feature.sub_features:
            if sub_feature.type != "CDS":
                continue
            chrom = record_id
            strand = sub_feature.strand
            CDS_start = sub_feature.location.start + 1 if strand == +1 else sub_feature.location.end
            CDS_end = sub_feature.location.end if strand == +1 else sub_feature.location.start + 1
            #region_start = CDS_start - (args.left * strand)
            #region_end = CDS_start + (args.right * strand)

            region_start_start = CDS_start - args.left if strand == +1 else CDS_start - args.right
            region_start_end = CDS_start + args.right if strand == +1 else CDS_start + args.left

            region_end_start = CDS_end - args.left if strand == +1 else CDS_end - args.right
            region_end_end = CDS_end + args.right if strand == +1 else CDS_end + args.left
            #print("aaa")
            start_coordinates = []
            end_coordinates = []
            for variant in variants:
                if record_id != variant.chrom:
                    continue
                if region_start_start <= variant.pos <= region_start_end:
                    start_coordinates.append((variant.pos - CDS_start) * strand)
                if region_end_start <= variant.pos <= region_end_end:
                    end_coordinates.append((variant.pos - CDS_end) * strand)
            all_variant_start_positions += start_coordinates
            all_variant_end_positions += end_coordinates
            #print(feature.qualifiers)
            gene_variants_positions.append([feature.qualifiers["Name"], strand, chrom, region_start_start,
                                            region_start_end, start_coordinates,
                                            region_end_start, region_end_end,
                                            end_coordinates])

with open(args.gene_variants, "w") as out_fd:
    out_fd.write("#gene\tstrand\tchrom\tregion_start_start\tregion_start_end\tstart_coordinates\tregion_end_start\tregion_end_end\tend_coordinates\n")
    for gene_annotations in gene_variants_positions:
        out_fd.write("\t".join([str(x) for x in gene_annotations[:-1]]) +
                     "\t%s\n" % (",".join([str(x) for x in gene_annotations[-1] ]) if gene_annotations[-1] else ".") )


bins = int((args.left + args.right) / 5)
plt.figure(1, dpi=150, figsize=(16, 16))

plt.subplot(2, 1, 1)
plt.hist(all_variant_start_positions, bins=bins)
plt.title("%i mutations" % len(all_variant_start_positions))
plt.xlim(xmin=-args.left, xmax=args.right)

plt.subplot(2, 1, 2)
plt.hist(all_variant_start_positions, bins=bins, normed=True)
plt.xlim(xmin=-args.left, xmax=args.right)

plt.savefig(args.start_histogramm_prefix + ".svg")
plt.savefig(args.start_histogramm_prefix + ".eps")
plt.close()

plt.figure(1, dpi=150, figsize=(16, 16))

plt.subplot(2, 1, 1)
plt.hist(all_variant_end_positions, bins=bins)
plt.title("%i mutations" % len(all_variant_end_positions))
plt.xlim(xmin=-args.left, xmax=args.right)

plt.subplot(2, 1, 2)
plt.hist(all_variant_end_positions, bins=bins, normed=True)
plt.xlim(xmin=-args.left, xmax=args.right)

plt.savefig(args.end_histogramm_prefix + ".svg")
plt.savefig(args.end_histogramm_prefix + ".eps")
plt.close()
