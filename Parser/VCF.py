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
from Parser.Abstract import Record, Collection
from Parser.CCF import RecordCCF, CollectionCCF


class RecordVCF(Record):
    def __init__(self, chrom, pos, id, ref, alt_list, qual, filter_list, info_dict, samples_list, description={}):
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
        #TODO: add data check
        #TODO: check parsing of files with several samples

    def __str__(self):
        #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample_1
        return self.string_form()

    def string_form(self):
        alt_string = ",".join(self.alt_list)
        filter_string = ";".join(self.filter_list)
        key_string_list = []
        for key in sorted(list(self.info_dict.keys())):
            if self.info_dict[key]:
                key_string_list.append(key + "=" + ",". join(map(lambda x: str(x), self.info_dict[key])))
            else:
                key_string_list.append(key)

        #[key + "=" + ",". join(map(lambda x: str(x), self.info_dict[key])) for key in sorted(list(self.info_dict.keys()))]
        info_string = ";".join(key_string_list)
        #format_string = ":".join(self.format_list)
        format_string = ":".join(self.samples_list[0].keys())
        samples_string = "\t".join([":".join([",".join(map(lambda x: str(x), sample[key])) for key in sample.keys()]) for sample in self.samples_list])

        #samples_string = "\t".join([":".join([",".join(str(sample[key])) for key in sample.keys()]) for sample in self.samples_list])

        return '\t'.join(map(lambda x: str(x), [self.chrom, self.pos, self.id, self.ref, alt_string,
                                                self.qual, filter_string, info_string, format_string, samples_string]))

    def check_indel(self):
        #checks if record is indel
        if len(self.ref) > 1 or len("".join(self.alt_list)) > len(self.alt_list):
            return True
        return False

    def gff_str(self, parent=None):
        # TODO: think how to rewrite, maybe remove
        attributes_string = "ID=Var%s%i" % (self.chrom, self.start)
        if parent:
            attributes_string += ";Parent=%s" % parent
        return "%s\tvariant_call\tvariant\t%i\t%i\t.\t.\t.\t%s" % (self.chrom, self.start, self.end, attributes_string)


class CollectionVCF(Collection):
    #TODO: rewrite metadata and header as classes to be consistent with abstract classes

    def __init__(self, metadata=None, record_list=None, header_list=None, vcf_file=None, from_file=True):
        self.linkage_dict = None
        if from_file:
            self.metadata = OrderedDict({})
            self.records = []
            with open(vcf_file, "r") as fd:
                for line in fd:
                    if line[:2] != "##":
                        self.header_list = line[1:].strip().split("\t")
                        self.samples = self.header_list[9:]
                        break
                    #print(line)
                    self._add_metadata(line)
                for line in fd:
                    self._add_record(line)
        else:
            self.metadata = metadata
            self.records = record_list
            self.header_list = header_list

    def __len__(self):
        return len(self.records)

    def _add_record(self, line):
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
            if self.metadata["INFO"][entry[0]]["Type"] == "Flag":
                info_dict[entry[0]] = []
            elif self.metadata["INFO"][entry[0]]["Type"] == "Integer":
                info_dict[entry[0]] = list(map(lambda x: int(x), entry[1].split(",")))
            elif self.metadata["INFO"][entry[0]]["Type"] == "Float":
                info_dict[entry[0]] = list(map(lambda x: float(x), entry[1].split(",")))
            else:
                info_dict[entry[0]] = entry[1].split(",")
        samples_list = []

        for sample_string in line_list[9:]:
            sample_dict = OrderedDict({})
            if sample_string == "./.":
                sample_dict["GT"] = "./."
                continue
            #else:
            for key, value_list in zip(line_list[8].split(":"), sample_string.split(":")):
                #print (key, value_list)
                if self.metadata["FORMAT"][key]["Type"] == "Integer":
                    sample_dict[key] = list(map(lambda x: int(x), value_list.split(",")))
                elif self.metadata["FORMAT"][key]["Type"] == "Float":
                    sample_dict[key] = list(map(lambda x: float(x), value_list.split(",")))
                else:
                    sample_dict[key] = value_list.split(",")
            #print(sample_dict)
            samples_list.append(sample_dict)
        #print(samples_list)
        self.records.append(RecordVCF(line_list[0], position, line_list[2], line_list[3],
                                      alt_list, quality, filter_list,
                                      info_dict, samples_list))

    def _split_by_equal_sign(self, string):
        try:
            index = string.index("=")
        except ValueError:
            #if "=" is not present in string (case of flag type in INFO field)
            return string, None
        return string[:index], string[index+1:]

    def _split_by_comma_sign(self, string):
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

    def _add_metadata(self, line):
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
            if key not in self.metadata:
                self.metadata[key] = OrderedDict({})
            self.metadata[key][value_id] = value

        else:
            self.metadata[key] = value

    def metadata2str(self):
        metadata_string = ""
        for key in self.metadata:
            if not isinstance(self.metadata[key], dict):
                metadata_string += "##%s=%s\n" % (key, self.metadata[key])
            else:
                prefix = "##%s=<" % key
                suffix = ">\n"
                for att_id in self.metadata[key]:
                    middle = "ID=%s," % att_id + ",".join(["%s=%s" % (param, self.metadata[key][att_id][param])
                                                           for param in self.metadata[key][att_id]])
                    metadata_string += prefix + middle + suffix
        return metadata_string

    def write(self, output_file):
        with open(output_file, "w") as out_fd:
            out_fd.write(self.metadata2str())
            out_fd.write("#" + "\t".join(self.header_list) + "\n")
            for record in self.records:
                out_fd.write(str(record) + "\n")

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
                                    header_list=self.header_list, from_file=False)
        heterozygotes = CollectionVCF(metadata=self.metadata, record_list=hetero_sites,
                                      header_list=self.header_list, from_file=False)
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

    def stack_regions(self):
        #Not implemented yet
        #TODO: write
        pass
        """
        #assume that records inside regions are sorted by start coordinate
        sequence_varcoord_dict = self.record_coordinates()
        staked = np.array([])
        if "sequence-region" in self.metadata:
            shift_dict = {}
            shift = 0
            for seqid in sorted(list(self.metadata["sequence-region"].keys())):
                shift_dict[seqid] = shift
                sequence_varcoord_dict[seqid] += shift - self.metadata["sequence-region"][seqid][0] + 1
                staked = np.hstack((staked, sequence_varcoord_dict[seqid]))
                length = self.metadata["sequence-region"][seqid][1] - \
                         self.metadata["sequence-region"][seqid][0] + 1
                shift += length
        else:
            #TODO:implement stacking for files without sequence-region metadata
            pass
        return {"All": staked}, shift_dict
        """

    def split_by_ref_and_alt(self, ref_alt_list):
        #TODO: check
        # structure of ref_alt_list:  [[ref1,[alt1.1, alt1.M1]], ..., [refN,[altN.1, ..., altN.MN]]]
        found_records = []
        filtered_out_records = []
        for record in self.records:
            #print (record.ALT)
            if (record.ref, record.alt_list) in ref_alt_list:
                found_records.append(record)
            else:
                filtered_out_records.append(record)
        found = CollectionVCF(metadata=self.metadata, record_list=found_records,
                                    header_list=self.header_list, from_file=False)
        filtered_out = CollectionVCF(metadata=self.metadata, record_list=filtered_out_records,
                                      header_list=self.header_list, from_file=False)
        return found, filtered_out

    def _split_ref(self, records):
        splited_dict = OrderedDict({"A": [], "C": [], "G": [], "T": [], "INDEL": []})
        nucleotides = ["A", "C", "G", "T"]
        for record in records:
            if record.ref in nucleotides:
                splited_dict[record.ref].append(record)
            else:
                splited_dict["INDEL"].append(record)
        return splited_dict

    def split_by_regions(self):
        #TODO: check
        regions_dict = self._split_regions()
        return [CollectionVCF(metadata=self.metadata, record_list=regions_dict[region],
                              header_list=self.header_list, from_file=False)
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
                          header_list=self.header_list, from_file=False).write(prefix + "_" + region + ".vcf")

    def get_positions(self):
        regions_dict = self._split_regions()
        positions_dict = {}
        for region in regions_dict:
            positions_dict[region] = np.array([record.pos for record in regions_dict[region]])
        return positions_dict

    def rainfall_plot(self, plot_name, base_colors=[], single_fig=True, dpi=150, figsize=(55, 70), facecolor="#D6D6D6",
                      ref_genome=None, min_gap_length=10, draw_gaps=False):
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
            fig.suptitle("Rainfall plot", fontsize=14)
            sub_plot_dict = OrderedDict({})
        index = 1
        for region in regions_dict:
            positions_dict[region] = np.array([record.pos for record in regions_dict[region]])
            #np.ediff1d return differences between consecutive elements in array, then 0 is added to the beginning
            distances_dict[region] = np.insert(np.ediff1d(positions_dict[region]),
                                               0, 0)
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
                for reference in region_reference_dict[region]:
                    plt.plot(region_reference_dict[region][reference][0],
                             region_reference_dict[region][reference][1],
                             color=reference_colors[reference],
                             marker='.', linestyle='None', label=reference)

                plt.title("Region %s" % region)
                #xlabel("Position")
                plt.ylabel("Distanse")
                plt.ylim(ymin=0)
                plt.axhline(y=100, color="#000000")
                plt.axhline(y=1000, color="#000000")
                plt.axhline(y=500, color="purple")
                plt.axhline(y=10, color="#000000")
                if ref_genome:
                    if draw_gaps:
                        for gap in ref_genome.gaps_dict[region]:
                            plt.gca().add_patch(plt.Rectangle((gap.location.start, 1),
                                                              gap.location.end - gap.location.start,
                                                              1024*32, facecolor="#aaaaaa", edgecolor='none'))
        if single_fig:

            plt.savefig("%s/%s.svg" % (plot_dir, plot_name))
            for region in sub_plot_dict:
                sub_plot_dict[region].set_yscale('log', basey=2)
                    #yscale ('log', basey = 2)
            plt.savefig("%s/%s_log_scale.svg" % (plot_dir, plot_name))
            plt.close()

    def hierarchical_clustering(self, method='average', dendrogramm_max_y=2000,
                                sample_name=None, save=False, clustering_dir="clustering",
                                dendrogramm_color_threshold=1000):
        # IMPORTANT! Use only for one-sample vcf
        # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage
        region_dict = self._split_regions()
        positions_dict = OrderedDict({})
        correlation_dict = OrderedDict({})
        linkage_dict = OrderedDict({})
        inconsistent_dict = OrderedDict({})
        clusters_dict = OrderedDict({})
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
            linkage_dict[region] = linkage(distance_matrix, method=method)
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
            np.savetxt("%s/inconsistent_coefficient_%s.t" % (clustering_dir, region), inconsistent_dict[region])

            #clusters_dict[region] = fcluster(linkage_dict[region], 1)
            #np.savetxt("clustering/clusters_%s.t" % region, clusters_dict[region], fmt="%i")
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
                     dendrogramm_color_threshold=1000):
        if self.linkage_dict:
            linkage_dict = self.linkage_dict
        else:
            region_dict, linkage_dict = self.hierarchical_clustering(method=cluster_distance,
                                                                     dendrogramm_max_y=dendrogramm_max_y,
                                                                     sample_name=sample_name,
                                                                     save=save_clustering,
                                                                     clustering_dir=clustering_dir,
                                                                     dendrogramm_color_threshold=dendrogramm_color_threshold)
        if split_by_regions:
            mut_clusters_dict = OrderedDict({})
        else:
            mut_clusters_list = []
        for region in linkage_dict:
            # http://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.fcluster.html#scipy.cluster.hierarchy.fcluster
            clusters_dict = OrderedDict({})
            clusters = fcluster(linkage_dict[region], threshold, criterion=extracting_method)
            for i in range(0, len(clusters)):
                if clusters[i] not in clusters_dict:
                    clusters_dict[clusters[i]] = [region_dict[region][i]]
                else:
                    clusters_dict[clusters[i]].append(region_dict[region][i])
            if split_by_regions:
                mut_clusters_dict[region] = \
                    CollectionCCF(record_list=[RecordCCF(vcf_records_list=clusters_dict[cluster],
                                                         from_records=True)
                                               for cluster in clusters_dict])
            else:
                mut_clusters_list += [RecordCCF(vcf_records_list=clusters_dict[cluster], from_records=True)
                                      for cluster in clusters_dict]
        if split_by_regions:
            return mut_clusters_dict
        return CollectionCCF(record_list=mut_clusters_list)

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
                for k in counted:
                    if k > 1:
                        nonsingleton += 1
                    if k > 2:
                        multicluster += 1
                n_nonsingleton_clusters.append(nonsingleton)
                n_multiclusters.append(multicluster)
            sub_plot_dict[region] = plt.subplot(side, side, index, axisbg="#D6D6D6")
            #ax = plt.gca()
            #ax.set_xticks(np.arange(0.5, 2.2, 0.1))

            plt.grid()
            plt.plot(coef_threshold_list, n_clusters_list, label="all")
            plt.plot(coef_threshold_list, n_nonsingleton_clusters, "green", label="2+")
            plt.plot(coef_threshold_list, n_multiclusters, "red", label="3+")
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
