#!/usr/bin/python2
import os
import re
from collections import OrderedDict
from math import sqrt

import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, inconsistent, cophenet, fcluster

import matplotlib.pyplot as plt

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation



#from General import check_path
from Parsers.Abstract import Record, Collection, Metadata, Header

#TODO: refactor whole file
ref_alt_variants = {"desaminases": [("C", ["T"]), ("G", ["A"])]
                    }


class RecordVCF(Record):
    def __init__(self, chrom, pos, id, ref, alt_list, qual, filter_list, info_dict, samples_list,
                 description={}, flags=None):
        self.chrom = chrom                              #str
        self.pos = pos                                  #int
        self.id = id                                    #str
        self.ref = ref                                  #str
        self.alt_list = alt_list                        #list, entries are strings
        self.qual = qual                                #real or "."
        self.filter_list = sorted(filter_list)          #list, entries are strings
        self.info_dict = info_dict                      #dict
        self.samples_list = samples_list                #list entries are dicts with keys from format_list and
                                                        #values are lists
        self.description = description
        self.flags = flags
        #TODO: add data check

    def __str__(self):
        alt_string = ",".join(self.alt_list)
        filter_string = ";".join(self.filter_list)
        key_string_list = []
        for key in sorted(list(self.info_dict.keys())):
            if self.info_dict[key]:
                key_string_list.append(key + "=" + ",". join(map(lambda x: str(x), self.info_dict[key])))
            else:
                key_string_list.append(key)

        info_string = ";".join(key_string_list)
        for sample in self.samples_list:
            if len(sample.keys()) > 1:
                format_string = ":".join(sample.keys())
                break

        samples_string = "\t".join([":".join([",".join(map(lambda x: str(x), sample[key])) for key in sample.keys()]) for sample in self.samples_list])
        return '\t'.join(map(lambda x: str(x), [self.chrom, self.pos, self.id, self.ref, alt_string,
                                                self.qual, filter_string, info_string, format_string, samples_string]))

    def check_indel(self):
        #checks if record is indel
        if len(self.ref) > 1 or len("".join(self.alt_list)) > len(self.alt_list):
            return True
        return False

    def count_samples(self):
        #counts samples with variants
        number = 0
        for sample in self.samples_list:
            if sample["GT"][0] != "./." and sample["GT"][0] != "0/0":
                number += 1
        return number

    def check_ref_alt_list(self, ref_alt_list, flag):
        if not self.flags:
            self.flags = set([])
        # structure of ref_alt_list:  [[ref1,[alt1.1, alt1.M1]], ..., [refN,[altN.1, ..., altN.MN]]]
        if (self.ref, self.alt_list) in ref_alt_list:
            self.flags.add(flag)

    def gff_str(self, parent=None):
        # TODO: think how to rewrite, maybe remove
        attributes_string = "ID=Var%s%i" % (self.chrom, self.start)
        if parent:
            attributes_string += ";Parent=%s" % parent
        return "%s\tvariant_call\tvariant\t%i\t%i\t.\t.\t.\t%s" % (self.chrom, self.start, self.end, attributes_string)

    def set_filter(self, expression, filter_name):
        if eval(expression):
            self.filter_list.append(filter_name)
            self.filter_list.sort()

    def add_info(self, info_name, info_value=None):
        value = info_value if isinstance(info_value, list) else [] if info_value is None else [info_value]
        if info_name in self.info_dict:
            self.info_dict[info_name] += value
        else:
            self.info_dict[info_name] = value

    def find_location(self, record_dict, key="Ftype", strand_key="Fstrand", genes_key="Genes",
                      genes_strand_key="Gstrand", feature_type_black_list=[],
                      use_synonym=False, synonym_dict=None, add_intergenic_label=True):
        # function is written for old variant (with sub_feature)s rather then new (with CompoundLocation)
        # id of one SeqRecord in record_dict must be equal to record.pos
        # locations will be written to description dictionary of record using "key" as key

        if key not in self.info_dict:
            self.info_dict[key] = set([])

        # strandness values:
        # N - not defined
        # B - both
        # P - plus
        # M - minus
        strands = ["B", "P", "M"]
        if strand_key not in self.info_dict:
            # by default
            self.info_dict[strand_key] = ["N"]
        for flag_key in (genes_key, genes_strand_key):
            if flag_key not in self.info_dict:
                self.info_dict[flag_key] = []

        #print(self.chrom, self.pos)
        for feature in record_dict[self.chrom].features:

            if feature.type in feature_type_black_list:
                continue

            if (self.pos - 1) in feature:
                if feature.type == "gene" or feature.type == "ncRNA":
                    #print(feature.qualifiers)
                    self.info_dict[genes_key].append(feature.qualifiers["Name"][0])
                    self.info_dict[genes_strand_key].append(strands[feature.strand])

                self.info_dict[key].add(self.get_synonym(feature.type, use_synonym=use_synonym,
                                                           synonym_dict=synonym_dict))
                if self.info_dict[strand_key][0] == "N":
                    self.info_dict[strand_key][0] = strands[feature.strand]
                elif strands[feature.strand] != self.info_dict[strand_key][0]:
                    self.info_dict[strand_key][0] = "B"
            else:
                continue

                #print(feature)
            for sub_feature in feature.sub_features:
                if sub_feature.type in feature_type_black_list:
                    continue
                if (self.pos - 1) in sub_feature:
                    self.info_dict[key].add(self.get_synonym(sub_feature.type, use_synonym=use_synonym,
                                                           synonym_dict=synonym_dict))
                    if self.info_dict[strand_key][0] == "N":
                        self.info_dict[strand_key][0] = strands[sub_feature.strand]
                    elif strands[sub_feature.strand] != self.info_dict[strand_key][0]:
                        self.info_dict[strand_key][0] = "B"

        if not self.info_dict[genes_key]:
            self.info_dict.pop(genes_key)
            self.info_dict.pop(genes_strand_key)

        if add_intergenic_label and (not self.info_dict[key]): # or ("gene" not in self.info_dict[key])):
            # igc == intergenic
            self.info_dict[key].add("igc")


class MetadataVCF(OrderedDict, Metadata):

    #def __init__(self, metadata=None):
    #    self.metadata = metadata if metadata is not None else OrderedDict({})

    @staticmethod
    def _split_by_equal_sign(string):
        try:
            index = string.index("=")
        except ValueError:
            #if "=" is not present in string (case of flag type in INFO field)
            return string, None
        return string[:index], string[index+1:]

    @staticmethod
    def _split_by_comma_sign(string):
        index_list = [-1]
        i = 1
        while (i < len(string)):
            if string[i] == "\"":
                i += 1
                while string[i] != "\"":
                    i += 1
            if string[i] == ",":
                index_list.append(i)
            i += 1
        index_list.append(len(string))
        return [string[index_list[j] + 1: index_list[j + 1]] for j in range(0, len(index_list) - 1)]

    def add_metadata(self, line):
        key, value = self._split_by_equal_sign(line[2:].strip())
        if value[0] == "<" and value[-1] == ">":
            #checking is value a list or no
            #print(value)
            #print(key, value)
            value = self._split_by_comma_sign(value[1:-1])
            #print(value)
            #parse in suppose that first parameter in value list is ID
            value_id = self._split_by_equal_sign(value[0])[1]
            #print(value[1:])
            value = OrderedDict(self._split_by_equal_sign(entry) for entry in value[1:])
            #print(value_id, value)
            if key not in self:
                self[key] = OrderedDict({})
            self[key][value_id] = value

        else:
            self[key] = value

    def __str__(self):
        metadata_string = ""
        for key in self:
            if not isinstance(self[key], dict):
                metadata_string += "##%s=%s\n" % (key, self[key])
            else:
                prefix = "##%s=<" % key
                suffix = ">\n"
                for att_id in self[key]:
                    middle = "ID=%s," % att_id + ",".join(["%s=%s" % (param, self[key][att_id][param])
                                                           for param in self[key][att_id]])
                    metadata_string += prefix + middle + suffix
        return metadata_string[:-1]


class HeaderVCF(list, Header):

    def __str__(self):
        return "#" + "\t".join(self)


class CollectionVCF(Collection):

    def __init__(self, metadata=None, record_list=None, header=None, vcf_file=None, samples=None, from_file=True, external_metadata=None):
        self.linkage_dict = None
        if from_file:
            self.metadata = MetadataVCF()
            self.records = []
            with open(vcf_file, "r") as fd:
                for line in fd:
                    if line[:2] != "##":
                        self.header = HeaderVCF(line[1:].strip().split("\t"))   #line[1:].strip().split("\t")
                        self.samples = self.header[9:]
                        break
                    #print(line)
                    self.metadata.add_metadata(line)
                for line in fd:
                    self.records.append(self.add_record(line, external_metadata=external_metadata))
        else:
            self.metadata = metadata
            self.records = [] if record_list is None else record_list
            self.header = header
            self.samples = samples

    def add_record(self, line, external_metadata=None):
        line_list = line.strip().split("\t")
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
        position = int(line_list[1])
        quality = "."
        if quality != line_list[5]:
            quality = float(line_list[5])
        alt_list = line_list[4].split(",")
        filter_list = line_list[6].split(",")          #list, entries are strings

        info_tuple_list = [self._split_by_equal_sign(entry) for entry in line_list[7].split(";")]
        info_dict = {}
        for entry in info_tuple_list:
            if self.metadata:
                if self.metadata["INFO"][entry[0]]["Type"] == "Flag":
                    info_dict[entry[0]] = []
                elif self.metadata["INFO"][entry[0]]["Type"] == "Integer":
                    info_dict[entry[0]] = list(map(lambda x: int(x), entry[1].split(",")))
                elif self.metadata["INFO"][entry[0]]["Type"] == "Float":
                    info_dict[entry[0]] = list(map(lambda x: float(x), entry[1].split(",")))
                else:
                    info_dict[entry[0]] = entry[1].split(",")
            else:
                if external_metadata["INFO"][entry[0]]["Type"] == "Flag":
                    info_dict[entry[0]] = []
                elif external_metadata["INFO"][entry[0]]["Type"] == "Integer":
                    info_dict[entry[0]] = list(map(lambda x: int(x), entry[1].split(",")))
                elif external_metadata["INFO"][entry[0]]["Type"] == "Float":
                    info_dict[entry[0]] = list(map(lambda x: float(x), entry[1].split(",")))
                else:
                    info_dict[entry[0]] = entry[1].split(",")
        samples_list = []

        for sample_string in line_list[9:]:
            #print (sample_string)
            sample_dict = OrderedDict({})
            if sample_string == "./.":
                sample_dict["GT"] = ["./."]
            else:
                for key, value_list in zip(line_list[8].split(":"), sample_string.split(":")):
                    #print (key, value_list)
                    if self.metadata:
                        if self.metadata["FORMAT"][key]["Type"] == "Integer":
                            sample_dict[key] = list(map(lambda x: int(x), value_list.split(",")))
                        elif self.metadata["FORMAT"][key]["Type"] == "Float":
                            sample_dict[key] = list(map(lambda x: float(x), value_list.split(",")))
                        else:
                            sample_dict[key] = value_list.split(",")
                    else:
                        if external_metadata["FORMAT"][key]["Type"] == "Integer":
                            sample_dict[key] = list(map(lambda x: int(x), value_list.split(",")))
                        elif external_metadata["FORMAT"][key]["Type"] == "Float":
                            sample_dict[key] = list(map(lambda x: float(x), value_list.split(",")))
                        else:
                            sample_dict[key] = value_list.split(",")
            #print(sample_dict)
            samples_list.append(sample_dict)
        #print(samples_list)
        return RecordVCF(line_list[0], position, line_list[2], line_list[3],
                                      alt_list, quality, filter_list,
                                      info_dict, samples_list)

    @staticmethod
    def _split_by_equal_sign(string):
        try:
            index = string.index("=")
        except ValueError:
            #if "=" is not present in string (case of flag type in INFO field)
            return string, None
        return string[:index], string[index+1:]

    def _split_by_comma_sign(self, string):
        return self._split_by_comma_sign(string, sign=",")

    @staticmethod
    def _split_by_sign(string, sign=","):
        # ignores sign in "
        index_list = [-1]
        i = 1
        while (i < len(string)):
            if string[i] == "\"":
                i += 1
                while string[i] != "\"":
                    i += 1
            if string[i] == sign:
                index_list.append(i)
            i += 1
        index_list.append(len(string))
        return [string[index_list[j] + 1: index_list[j + 1]] for j in range(0, len(index_list) - 1)]

    def split_by_zygoty(self):
        #splits sites by zygoty, site counts as heterozygote if even in one sample it it hetorozygote
        homo_sites = []
        hetero_sites = []
        for site in self.records:
            for sample_dict in site.samples_list:
                zyg = sample_dict["GT"][0].split("/")
                if zyg[0] != zyg[1]:
                    hetero_sites.append(site)
                    break
            else:
                homo_sites.append(site)
        homozygotes = CollectionVCF(metadata=self.metadata, record_list=homo_sites,
                                    header=self.header, from_file=False)
        heterozygotes = CollectionVCF(metadata=self.metadata, record_list=hetero_sites,
                                      header=self.header, from_file=False)
        return homozygotes, heterozygotes

    def record_coordinates(self, black_list=[], white_list=[]):
        #return dictionary, where keys are chromosomes and values numpy arrays of SNV coordinates
        sequence_varcoord_dict = {}
        for record in self.records:
            if black_list and (record.chrom in black_list):
                continue
            if white_list and (record.chrom not in white_list):
                continue
            if record.chrom not in sequence_varcoord_dict:
                sequence_varcoord_dict[record.chrom] = [record.pos]
            else:
                sequence_varcoord_dict[record.chrom].append(record.pos)
        for chrom in sequence_varcoord_dict:
            sequence_varcoord_dict[chrom] = np.array(sequence_varcoord_dict[chrom])
        return sequence_varcoord_dict

    def check_by_ref_and_alt(self, ref_alt_list, flag):
        for record in self:
            record.check_ref_alt_list(ref_alt_list, flag)

    def split_by_flags(self, flag_set, mode="all"):
        # possible modes:
        # all - record to be counted as 'with flag' must have all flags from flags_list
        # one - record to be counted as 'with flag' must have at least one flag from flags_list
        with_flag_records, without_flag_records = self.split_records_by_flags(flag_set, mode=mode)
        with_flag = CollectionVCF(metadata=self.metadata, record_list=with_flag_records,
                                  header=self.header, from_file=False)
        without_flag = CollectionVCF(metadata=self.metadata, record_list=without_flag_records,
                                     header=self.header, from_file=False)
        return with_flag, without_flag

    def split_by_ref_and_alt(self, ref_alt_list):
        # structure of ref_alt_list:  [[ref1,[alt1.1, alt1.M1]], ..., [refN,[altN.1, ..., altN.MN]]]
        found_records = []
        filtered_out_records = []
        for record in self.records:
            if (record.ref, record.alt_list) in ref_alt_list:
                found_records.append(record)
            else:
                filtered_out_records.append(record)
        found = CollectionVCF(metadata=self.metadata, record_list=found_records,
                              header=self.header, from_file=False)
        filtered_out = CollectionVCF(metadata=self.metadata, record_list=filtered_out_records,
                                     header=self.header, from_file=False)
        return found, filtered_out

    @staticmethod
    def _split_ref(records):
        splited_dict = OrderedDict({"A": [], "C": [], "G": [], "T": [], "INDEL": []})
        nucleotides = ["A", "C", "G", "T"]
        for record in records:
            if record.ref in nucleotides:
                splited_dict[record.ref].append(record)
            else:
                splited_dict["INDEL"].append(record)
        return splited_dict

    def set_filter(self, expression, filter_name):
        for record in self:
            if eval(expression):
                if "PASS" in record.filter_list or "." in record.filter_list:
                    record.filter_list = [filter_name]
                else:
                    record.filter_list.append(filter_name)

    def split_by_regions(self):
        #TODO: check
        regions_dict = self._split_regions()
        return [CollectionVCF(metadata=self.metadata, record_list=regions_dict[region],
                              header=self.header, from_file=False)
                for region in regions_dict]

    def _reference(self, record):
        nucleotides = ["A", "C", "G", "T"]
        if record.ref in nucleotides:
            return record.ref
        return "INDEL"

    def write_splited_regions(self, prefix):
        #TODO: check
        regions_dict = self._split_regions()

        for region in regions_dict:
            CollectionVCF(metadata=self.metadata, record_list=regions_dict[region],
                          header=self.header, from_file=False).write(prefix + "_" + region + ".vcf")

    def get_positions(self):
        regions_dict = self._split_regions()
        positions_dict = {}
        for region in regions_dict:
            positions_dict[region] = np.array([record.pos for record in regions_dict[region]])
        return positions_dict

    def rainfall_plot(self, plot_name, base_colors=[], single_fig=True, dpi=300, figsize=(40, 40), facecolor="#D6D6D6",
                      ref_genome=None, masked_regions=None, min_gap_length=10, draw_gaps=False, suptitle=None):
        print("Drawing rainfall plot...")
        plot_dir = "rainfall_plot"
        reference_colors = {"A": "#FBFD2B",    # yellow
                            "C": "#FF000F",     # red
                            "G": "#000FFF",     # blue
                            "T": "#4ED53F",     # green
                            "INDEL": "#000000"  # black
                            }
        if base_colors:
            reference_colors = base_colors


        regions_dict = self._split_regions()
        num_of_regions = len(regions_dict)
        positions_dict = {}
        distances_dict = {}
        region_reference_dict = {}
        os.system("mkdir -p %s" % plot_dir)
        if single_fig:
            fig = plt.figure(1, dpi=dpi, figsize=figsize, facecolor=facecolor)
            fig.suptitle(suptitle if suptitle else "Rainfall plot", fontsize=40, fontweight='bold', y=0.94)
            sub_plot_dict = OrderedDict({})
        index = 1
        for region in regions_dict:
            positions_dict[region] = np.array([record.pos for record in regions_dict[region]])
            #np.ediff1d return differences between consecutive elements in array, then 0 is added to the beginning
            distances_dict[region] = np.insert(np.ediff1d(positions_dict[region]), 0, 0)
            region_reference_dict[region] = OrderedDict({"A": [[], []],
                                                         "C": [[], []],
                                                         "G": [[], []],
                                                         "T": [[], []],
                                                         "INDEL": [[], []]})
            for i in range(0, len(regions_dict[region])):
                region_reference_dict[region][self._reference(regions_dict[region][i])][0].append(positions_dict[region][i])
                region_reference_dict[region][self._reference(regions_dict[region][i])][1].append(distances_dict[region][i])
            if single_fig:
                if not sub_plot_dict:
                    sub_plot_dict[region] = plt.subplot(num_of_regions, 1, index, axisbg=facecolor)
                else:
                    keys = list(sub_plot_dict.keys())
                    sub_plot_dict[region] = plt.subplot(num_of_regions, 1, index,
                                                        sharex=sub_plot_dict[keys[0]],
                                                        sharey=sub_plot_dict[keys[0]],
                                                        axisbg=facecolor)

                index += 1
                if draw_gaps:
                    if ref_genome:
                        for gap in ref_genome.gaps_dict[region]:
                            plt.gca().add_patch(plt.Rectangle((gap.location.start, 1),
                                                              gap.location.end - gap.location.start,
                                                              1024*32, facecolor="#777777", edgecolor='none'))
                    # masked regions should be SeqRecord dict
                    if masked_regions:
                        for feature in masked_regions[region].features:
                            plt.gca().add_patch(plt.Rectangle((int(feature.location.start)+1, 1),
                                                              feature.location.end - feature.location.start,
                                                              1024*32, facecolor="#aaaaaa", edgecolor='none'))

                for reference in region_reference_dict[region]:
                    plt.plot(region_reference_dict[region][reference][0],
                             region_reference_dict[region][reference][1],
                             color=reference_colors[reference],
                             marker='.', linestyle='None', label=reference)

                #plt.title("Region %s" % region)
                #xlabel("Position")
                plt.text(-0.08, 0.5, region, rotation=0, fontweight="bold", transform=sub_plot_dict[region].transAxes, fontsize=30,
                         horizontalalignment='center',
                         verticalalignment='center')
                plt.ylabel("Distanse")
                #plt.ylim(ymin=0)
                plt.axhline(y=100, color="#000000")
                plt.axhline(y=1000, color="#000000")
                plt.axhline(y=500, color="purple")
                plt.axhline(y=10, color="#000000")
                #if ref_genome:
                #    plt.xlim(xmax=len(ref_genome.reference_genome[region]))

        if single_fig:
            #plt.savefig("%s/%s.svg" % (plot_dir, plot_name))
            for region in sub_plot_dict:
                sub_plot_dict[region].set_yscale('log', basey=2)
                    #yscale ('log', basey = 2)
            plt.savefig("%s/%s_log_scale.svg" % (plot_dir, plot_name))
            plt.savefig("%s/%s_log_scale.pdf" % (plot_dir, plot_name))
            plt.savefig("%s/%s_log_scale.eps" % (plot_dir, plot_name))
            #plt.tight_layout()
            plt.close()

    def hierarchical_clustering(self, method='average', dendrogramm_max_y=2000,
                                sample_name=None, save=False, clustering_dir="clustering",
                                dendrogramm_color_threshold=1000,
                                draw_dendrogramm=True,
                                write_inconsistent=True,
                                write_correlation=True):
        # IMPORTANT! Use only for one-sample vcf
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
        region_dict = self._split_regions()
        positions_dict = OrderedDict({})
        correlation_dict = OrderedDict({})
        linkage_dict = OrderedDict({})
        inconsistent_dict = OrderedDict({})
        clusters_dict = OrderedDict({})
        if draw_dendrogramm or write_correlation or write_inconsistent:
            os.system("mkdir -p %s" % clustering_dir)
        for region in region_dict:
            positions_dict[region] = np.array([[record.pos] for record in region_dict[region]])

            # allowed methods(used to calculate distance between clusters):
            # 'complete'    -   Farthest Point Algorithm
            # 'single'      -   Nearest Point Algorithm
            # 'average'     -   UPGMA algorithm, distance between clusters is calculated as average from pairwise
            #                   distances between elements of clusters
            # 'weighted     -   WPGMA algorithm
            # 'centroid'    -   UPGMC algorithm
            # 'median'      -   WPGMC algorithm
            # 'ward'        -   incremental algorithm

            distance_matrix = pdist(positions_dict[region])
            #print(distance_matrix)
            linkage_dict[region] = linkage(distance_matrix, method=method)
            if draw_dendrogramm:
                plt.figure(1, dpi=150, figsize=(50, 20))
                dendrogram(linkage_dict[region],
                           color_threshold=dendrogramm_color_threshold,
                           leaf_font_size=4,
                           distance_sort=True)
                plt.ylim(ymax=dendrogramm_max_y)
                plt.axhline(y=500, color="purple")
                plt.axhline(y=1000, color="black")
                plt.savefig("%s/clustering_%s.svg" % (clustering_dir, region))
                plt.close()

            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.cophenet.html#scipy.cluster.hierarchy.cophenet
            # calculates cophenetic correlation coefficient to estimate accuracy of clustering
            correlation_dict[region] = cophenet(linkage_dict[region], distance_matrix)[0]

            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.inconsistent.html#scipy.cluster.hierarchy.inconsistent
            # calculates inconsistent coeff

            inconsistent_dict[region] = inconsistent(linkage_dict[region])
            if write_inconsistent:
                np.savetxt("%s/inconsistent_coefficient_%s.t" % (clustering_dir, region), inconsistent_dict[region])

            #clusters_dict[region] = fcluster(linkage_dict[region], 1)
            #np.savetxt("clustering/clusters_%s.t" % region, clusters_dict[region], fmt="%i")
        if write_correlation:
            sample = sample_name
            if not sample:
                sample = self.samples[0]
            with open("%s/correlation.t" % clustering_dir, "w") as cor_fd:
                cor_fd.write("sample\t%s\n" % ("\t".join(list(region_dict.keys()))))
                cor_fd.write("%s\t%s\n" % (sample, "\t".join([str(correlation_dict[region]) for region in region_dict])))

        if save:
            self.linkage_dict = linkage_dict

        return region_dict, linkage_dict

    def get_clusters(self,
                     extracting_method="inconsistent",
                     threshold=0.8,
                     cluster_distance='average',
                     dendrogramm_max_y=2000,
                     sample_name=None,
                     save_clustering=False,
                     clustering_dir="clustering",
                     split_by_regions=False,
                     dendrogramm_color_threshold=1000,
                     draw_dendrogramm=True,
                     return_collection=True,
                     write_inconsistent=True,
                     write_correlation=True):
        from Parsers.CCF import RecordCCF, CollectionCCF, MetadataCCF, HeaderCCF
        if self.linkage_dict:
            linkage_dict = self.linkage_dict
        else:
            region_dict, linkage_dict = self.hierarchical_clustering(method=cluster_distance,
                                                                     dendrogramm_max_y=dendrogramm_max_y,
                                                                     sample_name=sample_name,
                                                                     save=save_clustering,
                                                                     clustering_dir=clustering_dir,
                                                                     dendrogramm_color_threshold=dendrogramm_color_threshold,
                                                                     draw_dendrogramm=draw_dendrogramm,
                                                                     write_correlation=write_correlation,
                                                                     write_inconsistent=write_inconsistent)
        if split_by_regions:
            mut_clusters_dict = OrderedDict({})
        else:
            mut_clusters_list = []

        clusters = OrderedDict()
        for region in linkage_dict:
            clusters[region] = fcluster(linkage_dict[region], threshold, criterion=extracting_method)

        if return_collection:
            for region in region_dict:
                # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
                clusters_dict = OrderedDict({})

                for i in range(0, len(clusters[region])):
                    if clusters[region][i] not in clusters_dict:
                        clusters_dict[clusters[region][i]] = [region_dict[region][i]]
                    else:
                        clusters_dict[clusters[region][i]].append(region_dict[region][i])
                if split_by_regions:
                    mut_clusters_dict[region] = \
                        CollectionCCF(record_list=[RecordCCF(collection_vcf=CollectionVCF(record_list=clusters_dict[cluster], from_file=False),
                                                             from_records=True) for cluster in clusters_dict],
                                      metadata=MetadataCCF(self.samples, vcf_metadata=self.metadata, vcf_header=self.header),
                                      header=HeaderCCF("CLUSTER_ID\tCHROM\tSTART\tEND\tDESCRIPTION".split("\t")))
                else:
                    mut_clusters_list += [RecordCCF(collection_vcf=CollectionVCF(record_list=clusters_dict[cluster], from_file=False), from_records=True)
                                          for cluster in clusters_dict]
            if split_by_regions:
                return mut_clusters_dict
            return CollectionCCF(record_list=mut_clusters_list, metadata=MetadataCCF(self.samples, vcf_metadata=self.metadata, vcf_header=self.header),
                                 header=HeaderCCF("CLUSTER_ID\tCHROM\tSTART\tEND\tDESCRIPTION".split("\t")))
        else:
            return clusters

    def test_thresholds(self,
                        extracting_method="inconsistent",
                        threshold=None,
                        cluster_distance='average',
                        dendrogramm_max_y=2000,
                        sample_name=None,
                        save_clustering=False,
                        testing_dir="threshold_test"):

        # threshold is tuple(list) of three variables: min, max, number

        # extracting_method possible values
        #   inconsistent
        #   distance
        #   maxclust
        #   monocrit
        #   monocrit

        if self.linkage_dict:
            linkage_dict = self.linkage_dict
        else:
            region_dict, linkage_dict = self.hierarchical_clustering(method=cluster_distance,
                                                                     dendrogramm_max_y=dendrogramm_max_y,
                                                                     sample_name=sample_name,
                                                                     save=save_clustering,
                                                                     clustering_dir=testing_dir)

        num_of_regions = len(list(linkage_dict.keys()))

        side = int(sqrt(num_of_regions))
        if side*side != num_of_regions:
            side += 1
        sub_plot_dict = OrderedDict({})
        fig = plt.figure(2, dpi=150, figsize=(30, 30))
        fig.suptitle("Relashionship between number of clusters and threshold of %s" % extracting_method, fontsize=20)

        thresholds = threshold
        if extracting_method == "inconsistent":
            if not threshold:
                thresholds = (0.5, 1.5, 21)

        index = 1
        for region in linkage_dict:
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
            n_clusters_list = []
            n_nonsingleton_clusters = []
            n_multiclusters = []
            n_five_plus_clusters = []
            coef_threshold_list = np.linspace(*thresholds)  # best variant 0.5, 1.5, 21
            for i in coef_threshold_list:
                clusters = fcluster(linkage_dict[region], i, criterion=extracting_method)
                n_clusters_list.append(max(clusters))

                # counting clusters with 2+ and 3+ clusters
                ua, uind = np.unique(clusters, return_inverse=True)
                counted = np.bincount(uind)
                #j = 0
                nonsingleton = 0
                multicluster = 0  # 3+
                five_plus_clusters = 0 # 5+
                for k in counted:
                    if k > 1:
                        nonsingleton += 1
                    if k > 2:
                        multicluster += 1
                    if k > 4:
                        five_plus_clusters += 1
                n_nonsingleton_clusters.append(nonsingleton)
                n_multiclusters.append(multicluster)
                n_five_plus_clusters.append(five_plus_clusters)
            sub_plot_dict[region] = plt.subplot(side, side, index, axisbg="#D6D6D6")
            #ax = plt.gca()
            #ax.set_xticks(np.arange(0.5, 2.2, 0.1))

            plt.grid()
            plt.plot(coef_threshold_list, n_clusters_list, label="all")
            plt.plot(coef_threshold_list, n_nonsingleton_clusters, "green", label="2+")
            plt.plot(coef_threshold_list, n_multiclusters, "red", label="3+")
            plt.plot(coef_threshold_list, n_five_plus_clusters, "black", label="5+")
            plt.title("Region %s" % region)
            plt.legend(loc='upper right')
            plt.ylabel("Number of clusters")
            plt.xlabel("Threshold")
            #plt.axvline(x=0.8, color="purple")
            #plt.axvline(x=1.1, color="purple")

            plt.ylim(ymin=0)
            index += 1
        plt.savefig("%s/clusters_%s.svg" % (testing_dir, extracting_method))
        plt.close()

    def filter_by_expression(self, expression):
        filtered_records, filtered_out_records = self.filter_records_by_expression(expression)
        return CollectionVCF(metadata=self.metadata, record_list=filtered_records,
                               header=self.header, samples=self.samples, from_file=False), \
               CollectionVCF(metadata=self.metadata, record_list=filtered_out_records,
                               header=self.header, samples=self.samples, from_file=False)

    def add_info(self, metadata_line, expression, info_name, info_value=None):
        self.metadata.add_metadata(metadata_line)
        for record in self:
            if eval(expression):
                value = info_value if isinstance(info_value, list) else [] if info_value is None else [info_value]
                if info_name in record.info_dict:
                    record.info_dict[info_name] += value
                else:
                    record.info_dict[info_name] = value

    # methods for sum of two CollectionsVCF: no check for intersections(!!!!!!!!)
    def __add__(self, other):
        return CollectionVCF(metadata=self.metadata, record_list=self.records + other.records,
                             header=self.header, samples=self.samples, from_file=False)

    def __radd__(self, other):
        return CollectionVCF(metadata=self.metadata, record_list=self.records + other.records,
                             header=self.header, samples=self.samples, from_file=False)

    def parse_snpeff_info_record(self, string):
        effect, parameters = string.split("(")
        # remove closing bracket and split
        parameters = parameters[:-1].split("|")
        return [effect] + parameters

    def extract_snpeff_info(self, output_file):
        snpeff_info_dict_keys = "EFF", "LOS", "NMD"
        record_header_list = ["Chrom", "Pos", "Ref", "Alt"]
        snpeff_header_list = ["Effect", "Effect_Impact", "Functional_Class", "Codon_Change", "Amino_Acid_Change",
                  "Amino_Acid_Length", "Gene_Name", "Transcript_BioType", "Gene_Coding", "Transcript_ID",
                  "Exon_Rank", "Genotype_Number", "ERRORS", "WARNINGS"
                  ]

        #print(output_file)
        with open(output_file, "w") as out_fd:
            header_string = "#" + "\t".join(record_header_list + snpeff_header_list) + "\n"
            out_fd.write(header_string)
            for record in self:
                common_part = "%s\t%i\t%s\t%s" % (record.chrom, record.pos, record.ref, ",".join(record.alt_list))
                for effect in record.info_dict["EFF"]:
                    effect_parameters = self.parse_snpeff_info_record(effect)
                    num_parameters = len(effect_parameters)
                    for i in range(0, num_parameters):
                        if effect_parameters[i] == "":
                            effect_parameters[i] = "."
                    if num_parameters < 14:
                        effect_parameters += ["." for i in range(num_parameters, 14)]
                    out_fd.write(common_part + "\t" + "\t".join(effect_parameters) + "\n")

    def find_location(self, record_dict, key="Ftype", strand_key="Fstrand", genes_key="Genes", genes_strand_key="Gstrand",
                      feature_type_black_list=[],
                      use_synonym=False, synonym_dict=None, add_intergenic_label=True):

        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Types of features\">" % key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=1,Type=String,Description=\"Strand of features\">" % strand_key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Names of genes\">" % genes_key)
        self.metadata.add_metadata("##INFO=<ID=%s,Number=.,Type=String,Description=\"Strands of genes\">" % genes_strand_key)
        for record in self:
            record.find_location(record_dict, key=key, strand_key=strand_key,
                                 genes_key=genes_key, genes_strand_key=genes_strand_key,
                                 feature_type_black_list=feature_type_black_list,
                                 use_synonym=use_synonym, synonym_dict=synonym_dict,
                                 add_intergenic_label=add_intergenic_label)

    def count_strandness(self, prefix):
        count_dict = OrderedDict({})

        hor_coord_dict = {"C": 0, "G": 1}
        ver_coord_dict = {"N": 0, "P": 1, "M": 2, "B": 3}
        for record in self:
            if record.chrom not in count_dict:
                count_dict[record.chrom] = np.zeros((2, 4), dtype=int)
            count_dict[record.chrom][hor_coord_dict[record.ref]][ver_coord_dict[record.info_dict["Fstrand"][0]]] += 1

        count_dict["all"] = sum(count_dict.values())

        for chromosome in count_dict:
            with open("%s_%s.t" % (prefix, chromosome), "w") as out_fd:
                out_list = count_dict[chromosome].tolist()
                for index, name in zip(range(0, len(out_list)), ["C", "G"]):
                    out_list[index].insert(0, name)
                out_list.insert(0, [".", "N", "P", "M", "B"])
                for string_list in out_list:
                    out_fd.write("\t".join([str(x) for x in string_list]) + "\n")

        return count_dict

    def variants_start_end(self, left, right, record_dict, skip_genes_without_five_utr=False,
                           min_five_utr_len=10):

        gene_variants_positions = []
        all_variant_start_positions = []
        all_variant_end_positions = []

        for record_id in record_dict:
            for feature in record_dict[record_id].features:
                if feature.type != "gene":
                    continue
                if skip_genes_without_five_utr:
                    for sub_feature in feature.sub_features:
                        if sub_feature.type == "five_prime_UTR" and len(sub_feature) >= min_five_utr_len:
                            break
                    else:
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

                    region_start_start = CDS_start - left if strand == +1 else CDS_start - right
                    region_start_end = CDS_start + right if strand == +1 else CDS_start + left

                    region_end_start = CDS_end - left if strand == +1 else CDS_end - right
                    region_end_end = CDS_end + right if strand == +1 else CDS_end + left
                    #print("aaa")
                    start_coordinates = []
                    end_coordinates = []
                    for variant in self:
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
        return all_variant_start_positions, all_variant_end_positions, gene_variants_positions


class ReferenceGenome(object):
    # TODO: rewrite
    """docstring for ReferenceGenome"""
    chr_dict = {"chrI": "1",
                "chrII": "2",
                "chrIII": "3",
                "chrIV": "4",
                "chrV": "5",
                "chrVI": "6",
                "chrVII": "7",
                "chrVIII": "8",
                "chrIX": "9",
                "chrX": "10",
                "chrXI": "11",
                "chrXII": "12",
                "chrXIII": "13",
                "chrXIV": "14",
                "chrXV": "15",
                "chrXVI": "16",
                }
    feature_dict = {}
    gaps_dict = {}

    def __init__(self, ref_gen_file, index_file="refgen.idx", filetype="fasta"):
        self.ref_gen_file = ref_gen_file
        self.reference_genome = SeqIO.index_db(index_file, [ref_gen_file], filetype)
        """
        for entry in self.reference_genome:
            entry_name = self.chr_dict[self.reference_genome[entry].description.split(" ")[1]]
            self.feature_dict[entry_name] = []
            for feature in self.reference_genome[entry].features:
                if not (feature.type == "CDS" or feature.type == "source" or feature.type == "gene" or
                        feature.type == "rep_origin" or feature.type == "misc_feature"):
                    self.feature_dict[entry_name].append([feature.type,
                                                         [feature.location.start, feature.location.end,
                                                          feature.location.strand], feature.strand])
                #print feature.type,feature.location.start,feature.location.end, feature.location.strand
        """
    def find_gaps(self):
        gap_reg_exp = re.compile("N+", re.IGNORECASE)
        for region in self.reference_genome:
            self.gaps_dict[region] = []
            #print self.reference_genome[entry].seq
            gaps = gap_reg_exp.finditer(str(self.reference_genome[region].seq))  # iterator with
            for match in gaps:
                #print(match)
                self.gaps_dict[region].append(SeqFeature(FeatureLocation(match.start(), match.end()),
                                                         type="gap", strand=None))
        #print(self.gaps_dict)

if __name__ == "__main__":
    workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/LAN2xx/all"

    reference = ReferenceGenome("/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.fasta",
                                index_file="/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/LAN210_v0.9m.idx")
    #print(reference.reference_genome)
    reference.find_gaps()
    os.chdir(workdir)
    samples_list = sorted(os.listdir("."))

    suffix = "_GATK_best_merged.vcf"

    for sample in samples_list:
        print("Handling %s" % sample)

        os.chdir(workdir)
        os.chdir(sample)
        if "alignment_LAN210_v0.9m" not in os.listdir("."):
            continue
        os.chdir("alignment_LAN210_v0.9m")

        mutations = CollectionVCF(vcf_file=sample + suffix,
                                  from_file=True)
        print("Totaly %s mutations" % len(mutations))

        mutations.test_thresholds(save_clustering=True, testing_dir="testing_theshold_inconsistent")
        #mutations.test_thresholds(extracting_method='distance', threshold=(50, 5000, 100),
        #                          testing_dir="testing_threshold")
        #mutations.rainfall_plot(sample + suffix, ref_genome=reference, draw_gaps=True)
