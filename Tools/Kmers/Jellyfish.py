#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import numpy as np
from scipy.signal import argrelextrema

import matplotlib
matplotlib.use('Agg')
os.environ['MPLCONFIGDIR'] = '/tmp/'
import matplotlib.pyplot as plt
plt.ioff()
from Tools.Abstract import Tool
from Routines import MatplotlibRoutines, MathRoutines


class Jellyfish(Tool):
    """
    Several subcommands are not implemented: query, qhisto, qdump, qmerge, cite
    Not all options were implemented for count subcommand

    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "jellyfish", path=path, max_threads=max_threads)

    def count(self, in_file, out_file, kmer_length=23, hash_size=1000000, count_both_strands=False,
              lower_count=None, upper_count=None):
        # IMPORTANT! Not all options were implemented
        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = "-m %i" % kmer_length
        options += " -s %s" % hash_size
        options += " -t %i" % self.threads
        options += " -o %s" % out_file
        options += " -C" if count_both_strands else ""
        options += " -L %i" % lower_count if lower_count is not None else ""
        options += " -U %i" % upper_count if upper_count is not None else ""
        options += " %s" % (" ".join(in_file) if (isinstance(in_file, list)) or (isinstance(in_file, tuple)) else in_file)

        self.execute(options, cmd="jellyfish count")

    def stats(self, in_file, out_file, lower_count=None, upper_count=None):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = " -o %s" % out_file
        options += " -L %i" % lower_count if lower_count is not None else ""
        options += " -U %i" % upper_count if upper_count is not None else ""
        options += " %s" % in_file

        self.execute(options, cmd="jellyfish stats")

    def histo(self, in_file, out_file, bin_width=1, lower_count=1, upper_count=100000000,
              include_absent_kmers=False):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = " -o %s" % out_file
        options += " -t %i" % self.threads
        options += " -l %i" % lower_count
        options += " -h %i" % upper_count
        options += " -f" if include_absent_kmers else ""
        options += " -i %i" % bin_width
        options += " %s" % in_file

        self.execute(options, cmd="jellyfish histo")

    def dump(self, in_file, out_file, lower_count=None, upper_count=None,
             column_format=True, tab_separator=True):

        if (lower_count is not None) and (upper_count is not None):
            if lower_count > upper_count:
                raise ValueError("Upper limit for kmer counts is less than lower")

        options = " -o %s" % out_file
        options += " -L %i" % lower_count if lower_count is not None else ""
        options += " -U %i" % upper_count if upper_count is not None else ""
        options += " -c" if column_format else ""
        options += " -t" if tab_separator else ""
        options += " %s" % in_file

        self.execute(options, cmd="jellyfish dump")

    def get_kmer_list(self, in_file, out_prefix, kmer_length=23, hash_size=1000000, count_both_strands=False,
                      lower_count=None, upper_count=None):
        base_file = "%s_%i_mer.jf" % (out_prefix, kmer_length)
        kmer_table_file = "%s_%i_mer.counts" % (out_prefix, kmer_length)
        kmer_file = "%s_%i_mer.kmer" % (out_prefix, kmer_length)
        self.count(in_file, base_file, kmer_length=kmer_length, hash_size=hash_size,
                   count_both_strands=count_both_strands,
                   lower_count=lower_count, upper_count=upper_count)
        self.dump(base_file, kmer_table_file)
        sed_string = 'sed -e "s/\t.*//" %s > %s' % (kmer_table_file, kmer_file)
        os.system(sed_string)

    def draw_kmer_distribution(self, histo_file, kmer_length, output_prefix, output_formats=["svg", "png"],
                               logbase=10, non_log_low_limit=5, non_log_high_limit=100, order=3, mode="wrap",
                               check_peaks_coef=10, draw_separated_pictures=False):
        bins, counts = np.loadtxt(histo_file, unpack=True)
        maximums_to_show, minimums_to_show, \
            unique_peak_borders, number_of_distinct_kmers, \
            number_of_distinct_kmers_with_errors,\
            total_number_of_kmers, total_number_of_kmers_with_errors, \
            estimated_genome_size = self.extract_parameters_from_histo(counts, bins,
                                                                             output_prefix,
                                                                             order=order,
                                                                             mode=mode,
                                                                             check_peaks_coef=check_peaks_coef)

        unique_peak_width = unique_peak_borders[1] - unique_peak_borders[0] + 1
        #print unique_peak_borders
        #print unique_peak_width
        #print "Maximums to show"
        #print maximums_to_show
        #print "Minimums to show"
        #print minimums_to_show

        unique_peak_borders_mean_multiplicity = MathRoutines.mean_from_bins(bins[unique_peak_borders[0]: unique_peak_borders[1]+1],
                                                                            counts[unique_peak_borders[0]: unique_peak_borders[1]+1])
        std_1 = MathRoutines.std_from_bins(bins[unique_peak_borders[0]: unique_peak_borders[1]+1],
                                           counts[unique_peak_borders[0]: unique_peak_borders[1]+1],
                                           mean=unique_peak_borders_mean_multiplicity)
        var_1 = std_1 / unique_peak_borders_mean_multiplicity

        fraction_of_distinct_kmers_with_errors = float(number_of_distinct_kmers_with_errors)/float(number_of_distinct_kmers)
        fraction_of_kmers_with_errors = float(total_number_of_kmers_with_errors)/float(total_number_of_kmers)
        general_stats = "Number of distinct kmers\t%i\n" % number_of_distinct_kmers
        general_stats += "Number of distinct kmers with errors\t%i\n" % number_of_distinct_kmers_with_errors
        general_stats += "Fraction of distinct kmers with errors\t%.3f\n" % np.around(fraction_of_distinct_kmers_with_errors, decimals=3)

        general_stats += "Total number of kmers\t%i\n" % total_number_of_kmers
        general_stats += "Total number of kmers with errors\t%i\n" % total_number_of_kmers_with_errors
        general_stats += "Fraction of kmers with errors\t%.3f\n" % np.around(fraction_of_kmers_with_errors, decimals=3)
        general_stats += "Kmer multiplicity at first minimum\t%s\n" % (str(minimums_to_show[0][0]) if minimums_to_show else "None")
        general_stats += "Kmer multiplicity at first maximum\t%s\n" % (str(maximums_to_show[0][0]) if maximums_to_show else "None")
        general_stats += "Width of first peak\t%i\n" % unique_peak_width
        general_stats += "Mean kmer multiplicity in first peak\t%.2f\n" % np.around(unique_peak_borders_mean_multiplicity, decimals=2)
        general_stats += "Standard deviation of kmer multiplicity in first peak\t%.2f\n" % np.around(std_1, decimals=2)
        general_stats += "Variance coefficient of kmer multiplicity in first peak\t%.2f\n" % np.around(var_1,
                                                                                                       decimals=2)
        general_stats += "Estimated genome size, bp\t%i\n" % estimated_genome_size
        with open("%s.histo.stats" % output_prefix, "w") as stat_fd:
            stat_fd.write(general_stats)
        print(general_stats)

        max_bin = max(bins)

        selected_counts = counts[non_log_low_limit-1:non_log_high_limit]
        selected_bins = bins[non_log_low_limit-1:non_log_high_limit]

        if draw_separated_pictures:
            figure = plt.figure(1, figsize=(8, 8), dpi=300)
            subplot = plt.subplot(1, 1, 1)
            plt.suptitle("Distribution of %i-mers" % kmer_length, fontweight='bold')
            plt.plot(bins, counts)
            plt.xlim(xmin=1, xmax=max_bin)
            plt.xlabel("Multiplicity")
            plt.ylabel("Number of distinct %s-mers" % kmer_length)
            subplot.set_yscale('log', basey=logbase)
            subplot.set_xscale('log', basex=logbase)

            for extension in output_formats:
                plt.savefig("%s.logscale.%s" % (output_prefix, extension))

            plt.close()

            figure = plt.figure(2, figsize=(8, 8), dpi=300)
            subplot = plt.subplot(1, 1, 1)
            plt.suptitle("Distribution of %s-mers" % kmer_length, fontweight='bold')
            plt.plot(selected_bins, selected_counts)

            plt.xlabel("Multiplicity")
            plt.ylabel("Number of distinct %s-mers" % kmer_length)
            plt.xlim(xmin=non_log_low_limit, xmax=non_log_high_limit)

            for extension in output_formats:
                plt.savefig("%s.no_logscale.%s" % (output_prefix, extension))

            plt.close()

        size_in_gigabases = float(estimated_genome_size) / float(10 ** 9)
        size_in_megabases = float(estimated_genome_size) / float(10 ** 6)
        if size_in_gigabases > 0:
            legend = "Genome size %.2f G" % size_in_gigabases
        else:
            legend = "Genome size %.2f M" % size_in_megabases

        for index in range(3, 5):
            figure = plt.figure(index, figsize=(6, 12), dpi=400)
            subplot_list = []
            for i, b, c in zip([1, 2], [bins, selected_bins], [counts, selected_counts]):
                subplot_list.append(plt.subplot(2, 1, i))
                plt.suptitle("Distribution of %s-mers" % kmer_length, fontweight='bold', fontsize=13)
                plt.plot(b, c)

                plt.legend((legend, ), loc="upper right")

                if index == 4:
                    for minimum in minimums_to_show:
                        plt.plot([minimum[0], minimum[0]], [0, minimum[1]], 'r--', lw=1)
                    for maximum in maximums_to_show:
                        plt.plot([maximum[0], maximum[0]], [0, maximum[1]], 'g--', lw=1)

                plt.ylabel("Number of distinct %s-mers" % kmer_length, fontsize=13)
                if i == 1:
                    subplot_list[0].set_yscale('log', basey=logbase)
                    subplot_list[0].set_xscale('log', basex=logbase)
                    plt.xlim(xmin=1, xmax=max_bin)
                elif i == 2:
                    plt.xlim(xmin=non_log_low_limit, xmax=non_log_high_limit)
                    plt.xlabel("Multiplicity", fontsize=15)

            MatplotlibRoutines.zoom_effect(subplot_list[0], subplot_list[1], non_log_low_limit, non_log_high_limit)
            plt.subplots_adjust(hspace=0.12, wspace=0.05, top=0.95, bottom=0.05, left=0.14, right=0.95)

            for extension in output_formats:
                plt.savefig("%s.%s" % ("%s.peaks_and_gaps" % output_prefix if index == 4 else output_prefix, extension))

            plt.close()

    @staticmethod
    def find_peak_indexes_from_histo(counts, order=3, mode="wrap"):
        """
        order:
            How many points on each side to use for the comparison to consider comparator(n, n+x) to be True.
        mode:
            How the edges of the vector are treated. 'wrap' (wrap around) or 'clip' (treat overflow as the same
            as the last
        """
        local_maximums_idx = argrelextrema(counts, np.greater, order=order, mode=mode)[0]
        local_minimums_idx = argrelextrema(counts, np.less, order=order, mode=mode)[0]

        return local_minimums_idx, local_maximums_idx

    @staticmethod
    def extract_parameters_from_histo(counts, bins, output_prefix, order=3, mode="wrap", check_peaks_coef=10):
        """
        check_peaks_coef:
            histogram is checked for presence of additional peaks in range [first_unique_peak, check_peaks_coef*first_unique_peak]
        """

        local_maximums_idx = argrelextrema(counts, np.greater, order=order, mode=mode)[0]
        local_minimums_idx = argrelextrema(counts, np.less, order=order, mode=mode)[0]

        with open("%s.local_maximums" % output_prefix, "w") as out_fd:
            out_fd.write("#multiplicity\tnumber_of_kmers\n")
            for idx in local_maximums_idx:
                out_fd.write("%i\t%i\n" % (bins[idx], counts[idx]))

        with open("%s.local_minimums" % output_prefix, "w") as out_fd:
            out_fd.write("#multiplicity\tnumber_of_kmers\n")
            for idx in local_minimums_idx:
                out_fd.write("%i\t%i\n" % (bins[idx], counts[idx]))

        first_unique_peak_idx_idx = 0 if local_maximums_idx[0] != 0 else 1
        first_unique_peak_coverage = bins[local_maximums_idx[first_unique_peak_idx_idx]]

        max_checked_coverage = check_peaks_coef * first_unique_peak_coverage
        peaks_in_checked_area_idx = [local_maximums_idx[first_unique_peak_idx_idx]]
        minimums_in_checked_area_idx = []

        for i in range(first_unique_peak_idx_idx+1, len(local_maximums_idx)):
            if bins[local_maximums_idx[i]] <= max_checked_coverage:
                peaks_in_checked_area_idx.append(local_maximums_idx[i])
            else:
                break

        for minimum_index in local_minimums_idx:
            if bins[minimum_index] <= max_checked_coverage:
                minimums_in_checked_area_idx.append(minimum_index)

        if len(peaks_in_checked_area_idx) > 1:
            print("WARNING! Additional k-mer peaks were detected with multiplicity (%i, %i]" % (first_unique_peak_coverage,
                                                                                                max_checked_coverage))

        nearest_value_to_first_min_idx = MathRoutines.find_nearest_scalar(counts[local_maximums_idx[first_unique_peak_idx_idx]:],
                                                                          counts[local_minimums_idx[0]]) + local_maximums_idx[first_unique_peak_idx_idx]

        number_of_distinct_kmers = sum(counts)
        number_of_distinct_kmers_with_errors = sum(counts[0:local_minimums_idx[0]])
        total_number_of_kmers = sum(np.multiply(counts, bins))
        total_number_of_kmers_with_errors = sum(np.multiply(counts[0:local_minimums_idx[0]],
                                                                     bins[0:local_minimums_idx[0]]))

        maximums_to_show = [(bins[i], counts[i]) for i in peaks_in_checked_area_idx]
        minimums_to_show = [(bins[i], counts[i]) for i in minimums_in_checked_area_idx]

        estimated_genome_size = 0
        if minimums_to_show:
            for i in range(int(minimums_to_show[0][0]), len(counts)):
                estimated_genome_size += counts[i] * bins[i]
        else:
            print("WARNING: Minimum between error peak and unique peak was not found. "
                  "Estimation of genome size is impossible ")
            estimated_genome_size = None

        estimated_genome_size = estimated_genome_size/first_unique_peak_coverage if estimated_genome_size else None

        return maximums_to_show, \
               minimums_to_show, \
               (local_minimums_idx[0], nearest_value_to_first_min_idx), \
               number_of_distinct_kmers, number_of_distinct_kmers_with_errors, \
               total_number_of_kmers, total_number_of_kmers_with_errors, estimated_genome_size

