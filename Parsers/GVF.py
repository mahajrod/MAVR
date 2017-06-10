__author__ = 'mahajrod'

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import matplotlib.ticker

from matplotlib.patches import Rectangle
from math import log

#TODO: rewrite for less memomory usage


class RecordGVF():
    def __init__(self, seqid, source, type, start, end, score, strand, phase, attributes_dict):
        self.seqid = seqid          #str
        self.source = source        #str
        self.type = type            #str
        self.start = start          #int
        self.end = end              #int
        self.score = score          #real or "."
        self.strand = strand        #"+" or "-"
        self.phase = phase          #str
        if ("ID" not in attributes_dict) or ("Variant_seq" not in attributes_dict) or ("Reference_seq" not in attributes_dict):
            raise ValueError("Attribute dictionary doesnt have one or more required keys ('ID', 'Variant_seq', 'Reference_seq')")
        self.attributes_dict = attributes_dict #dict

    def __str__(self):
        return self.string_form()

    def string_form(self):
        attributes_str = ';'.join([key + "=" + self.attributes_dict[key] for key in self.attributes_dict])
        return '\t'.join(map(lambda x: str(x), [self.seqid, self.source, self.type, self.start, self.end,
                                                self.score, self.strand, self.phase, attributes_str]))


class CollectionGVF():
    def __add_record(self, line):
        line_list = line.strip().split("\t")
        line_list[3] = int(line_list[3])
        line_list[4] = int(line_list[4])
        if line_list[5] != ".":
            line_list[5] = float(line_list[5])
        line_list[8] = {key: value for key, value in map(lambda x: x.split("="), line_list[8].split(";"))}
        self.records.append(RecordGVF(*line_list))

    def __add_metadata(self, line):
        tmp = line[2:].strip().split(" ")
        key, value = tmp[0], " ".join(tmp[1:])
        if key in self.metadata:
            self.metadata[key].append(value)
        else:
            self.metadata[key] = [value]

    def __init__(self, metadata=None, record_list=None, gvf_file=None, from_file=True):
        if not from_file:
            self.metadata = metadata
            self.records = record_list
        else:
            self.metadata = {}
            self.records = []
            with open(gvf_file, "r") as in_fd:
                for line in in_fd:
                    if line[0] != "#":
                        self.__add_record(line)
                        break
                    self.__add_metadata(line)
                for line in in_fd:
                    self.__add_record(line)
        self.longest_region_length = None
        if "sequence-region" in self.metadata:
            tmp_list = []
            for i in range(0, len(self.metadata["sequence-region"])):
                tmp = self.metadata["sequence-region"][i].split(" ")
                tmp_list.append((tmp[0], list(map(lambda x: int(x), tmp[1:]))))
            self.metadata["sequence-region"] = {key: value for key, value in tmp_list}

            self.longest_region_length = max([self.metadata["sequence-region"][region][1] -
                                              self.metadata["sequence-region"][region][0] + 1
                                              for region in self.metadata["sequence-region"]])

    def record_coordinates(self):
        #return dictionary, where keys are seqids and values numpy arrays of SNV coordinates
        sequence_varcoord_dict = {}
        for record in self.records:
            if record.seqid not in sequence_varcoord_dict:
                sequence_varcoord_dict[record.seqid] = [record.start]
            else:
                sequence_varcoord_dict[record.seqid].append(record.start)
        for seqid in sequence_varcoord_dict:
            sequence_varcoord_dict[seqid] = np.array(sequence_varcoord_dict[seqid])
        return sequence_varcoord_dict

    def stack_regions(self):
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

    def rainfall_plot(self, output_file):
        #TODO: write rainfall plot
        pass

    def coverage_plot(self, cov_outfile="variation_coverage_plot.svg", distr_outfile="variation_distobution_plot.svg",
                      subplots=None,
                      cov_figsize=(10, 10), cov_dpi =300, cov_nbins=1000, cov_ylabel="N of SNV",
                      cov_xlabel="Position", share_x=False, share_y=True, single_plot=False,
                      distr_figsize=(10, 10), distr_dpi =300, distr_nbins=None, distr_ylabel="N of windows",
                      distr_xlabel="N of SNV in window",
                      black_list=[], white_list=[]):
        shift_dict = {}
        if single_plot:
            sequence_varcoord_dict, shift_dict = self.stack_regions()
        else:
            sequence_varcoord_dict = self.record_coordinates()
        num_of_seq = len(sequence_varcoord_dict)

        if (not subplots) or single_plot:
            subplot_tuple = (num_of_seq, 1)
        else:
            subplot_tuple = subplots

        #plt.figure(1, figsize=figsize, dpi=dpi)
        cov_fig, cov_axes = plt.subplots(nrows=subplot_tuple[0], ncols=subplot_tuple[1],
                                         figsize=cov_figsize, dpi=cov_dpi)
        #fig.suptitle('SNV coverage', fontsize=14, fontweight='bold')
        if not single_plot:
            cov_fig.tight_layout()
        else:
            cov_fig.suptitle('SNV coverage', fontweight='bold')
        index = 1
        subplot_list = []
        hist_dict = {}          #dict of n,bins,patches for each histogramm
        for seqid in sorted(list(sequence_varcoord_dict.keys())):
            if black_list and (seqid in black_list):
                continue
            if white_list and (seqid not in white_list):
                continue

            if subplot_list:
                subplot_list.append(plt.subplot(subplot_tuple[0], subplot_tuple[1], index, axisbg="#aaaaaa",
                                    sharey=subplot_list[0]))#, sharex=subplot_list[0]))
            else:
                subplot_list.append(plt.subplot(subplot_tuple[0], subplot_tuple[1], index, axisbg="#aaaaaa"))

            n_of_bins = cov_nbins
            if not single_plot:
                if "sequence-region" in self.metadata:
                    #plt.gca().add_patch(Rectangle((0, 0), self.metadata["sequence-region"][seqid][1], 900,facecolor="#aaaaaa", edgecolor='none'))
                    subplot_list[index-1].set_xlim(self.metadata["sequence-region"][seqid][0],
                                                   self.metadata["sequence-region"][seqid][1])

                if self.longest_region_length:
                    n_of_bins = 1 + int(cov_nbins * self.metadata["sequence-region"][seqid][1]/self.longest_region_length)
            else:
                last_region = sorted(list(self.metadata["sequence-region"].keys()))[-1]
                subplot_list[index-1].set_xlim(0, shift_dict[last_region] +
                                                  self.metadata["sequence-region"][last_region][1])

            hist_dict[seqid] = plt.hist(sequence_varcoord_dict[seqid], bins=n_of_bins, color="green")

            formatter = matplotlib.ticker.ScalarFormatter()
            formatter.set_powerlimits((-3, 2))
            subplot_list[index-1].xaxis.set_major_formatter(formatter)
            subplot_list[index-1].yaxis.set_major_formatter(formatter)
            plt.title("Chromosome %s" % seqid)
            #plt.legend(loc='upper right')
            plt.ylabel(cov_ylabel)
            #plt.xlabel(cov_xlabel)
            index += 1
        if shift_dict:
            #add chromosome edges
            for seqid in shift_dict:
                plt.plot((shift_dict[seqid], shift_dict[seqid]), (0, 400), 'k-', color="red")
            #plt.legend(loc=
        #print(hist_dict)
        plt.savefig(cov_outfile)
        plt.close()

        distr_fig, distr_axes = plt.subplots(nrows=subplot_tuple[0], ncols=subplot_tuple[1],
                                             figsize=distr_figsize, dpi=distr_dpi)
        if not single_plot:
            distr_fig.tight_layout()
        else:
            cov_fig.suptitle('SNV coverage distrobution', fontweight='bold')
        index = 1
        distr_subplot_list = []
        for seqid in sorted(list(hist_dict.keys())):
            """
            if black_list and (seqid in black_list):
                continue
            if white_list and (seqid not in white_list):
                continue
            """
            if distr_subplot_list:
                distr_subplot_list.append(plt.subplot(subplot_tuple[0], subplot_tuple[1], index, axisbg="#aaaaaa",
                                          sharey=subplot_list[0]))#, sharex=subplot_list[0]))
            else:
                distr_subplot_list.append(plt.subplot(subplot_tuple[0], subplot_tuple[1], index, axisbg="#aaaaaa"))

            n_of_bins = distr_nbins
            if not n_of_bins:
                n_of_bins = 1 + 3.32 * log(len(hist_dict[seqid][0]), 10)

            plt.hist(hist_dict[seqid][0], bins=n_of_bins, color="green")

            formatter = matplotlib.ticker.ScalarFormatter()
            formatter.set_powerlimits((-3, 2))
            subplot_list[index-1].xaxis.set_major_formatter(formatter)
            subplot_list[index-1].yaxis.set_major_formatter(formatter)
            plt.title("Chromosome %s" % seqid)
            #plt.legend(loc='upper right')
            plt.ylabel(distr_ylabel)
            plt.xlabel(distr_xlabel)
            #plt.xlabel(cov_xlabel)
            index += 1
        plt.savefig(distr_outfile)
        plt.close()

if __name__ == "__main__":
    gvf_col = CollectionGVF(gvf_file="/home/mahajrod/genetics/tetraodon/SNV/Tetraodon_nigroviridis.gvf", from_file=True)
    #print(gvf_col.metadata["sequence-region"])
    #gvf_col.coverage_plot(figsize=(30, 15), nbins=3000, subplots=(6, 5))#black_list=["Un_random", "1_random",
                                                                         #"15_random", "2_random",
                                                                         #"21_random"])
    gvf_col.coverage_plot(cov_outfile="variation_coverage_plot_single.svg",
                          distr_outfile="variation_distrobution_plot_single.svg",
                          cov_figsize=(22.5, 7.5), cov_nbins=3000,
                          distr_figsize=(10, 10), distr_nbins=50, single_plot=True)







