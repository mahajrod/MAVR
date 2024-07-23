#!/usr/bin/env python
"""
TADbit based code. It contains slightly modified quality_plot function from TADbit.
It was modified to produce machine-readable file with stats including percentage of diggested sites.
Also, a generation of a figure in multiple formats was added as option.
"""

import argparse

from warnings                             import catch_warnings, simplefilter
from itertools                            import zip_longest
import re

from collections                          import OrderedDict
from numpy                                import nanstd, nanmean, linspace, nansum

import matplotlib.pyplot as plt

from pytadbit.utils.file_handling         import magic_open
from pytadbit.mapping.restriction_enzymes import RESTRICTION_ENZYMES, iupac2regex
from pytadbit.mapping.restriction_enzymes import religateds, repaired

try:
    basestring
except NameError:
    basestring = str


def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def quality_plot(fnam, r_enz=None, nreads=float('inf'), axe=None,
                 output_prefix=None, extension_list=["png", "svg"], paired=False):
    """
    Plots the sequencing quality of a given FASTQ file. If a restrinction enzyme
    (RE) name is provided, can also represent the distribution of digested and
    undigested RE sites and estimate an expected proportion of dangling-ends.

    Proportion of dangling-ends is inferred by counting the number of times a
    dangling-end site, is found at the beginning of any of the reads (divided by
    the number of reads).

    :param fnam: path to FASTQ file
    :param None nreads: max number of reads to read, not necesary to read all
    :param None savefig: path to a file where to save the image generated;
       if None, the image will be shown using matplotlib GUI (the extension
       of the file name will determine the desired format).
    :param False paired: is input FASTQ contains both ends

    :returns: the percentage of dangling-ends (sensu stricto) and the percentage of
       reads with at least a ligation site.
    """
    phred = dict([(c, i) for i, c in enumerate(
        '!"#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~')])
    if isinstance(r_enz, list):
        r_enzs = r_enz
    elif isinstance(r_enz, basestring):
        r_enzs = [r_enz]
    for k in list(RESTRICTION_ENZYMES.keys()):
        for i in range(len(r_enzs)):
            if k.lower() == str(r_enz[i]).lower():
                r_enz[i] = k
    # else let it as None

    quals = []
    henes = []
    sites = {}
    fixes = {}
    liges = OrderedDict()
    ligep = {}
    tkw = dict(size=4, width=1.5)
    fhandler = magic_open(fnam)
    if len(r_enzs) == 1 and r_enzs[0] is None:
        if nreads:
            while True:
                try:
                    next(fhandler)
                except EOFError:
                    break
                seq = next(fhandler)
                if 'N' in seq:
                    henes.extend([i for i, s in enumerate(seq) if s == 'N'])
                next(fhandler)
                line = next(fhandler)
                quals.append([phred[i] for i in line.strip()])
                if len(quals) > nreads:
                    break
        else: # do this because it's faster
            while True:
                try:
                    next(fhandler)
                except EOFError:
                    break
                seq = next(fhandler)
                if 'N' in seq:
                    henes.extend([i for i, s in enumerate(seq) if s == 'N'])
                next(fhandler)
                line = next(fhandler)
                quals.append([phred[i] for i in line.strip()])
    else:
        r_sites = {}
        d_sites = {}
        for r_enz in r_enzs:
            r_sites[r_enz] = RESTRICTION_ENZYMES[r_enz].replace('|', '')
            d_sites[r_enz] = repaired(r_enz)
            sites[r_enz] = []  # initialize dico to store undigested sites
            fixes[r_enz] = []  # initialize dico to store digested sites
        l_sites = religateds(r_enzs)
        l_sites = OrderedDict((k,iupac2regex(l_sites[k])) for k in l_sites)
        site = {}
        fixe = {}
        for r_enz in r_enzs:
            site[r_enz] = re.compile(iupac2regex(r_sites[r_enz]))
            fixe[r_enz] = re.compile(iupac2regex(d_sites[r_enz]))
        # ligation sites should appear in lower case in the sequence
        lige = {}
        for k in l_sites:
            liges[k] = []  # initialize dico to store sites
            ligep[k] = 0   # initialize dico to store sites
            l_sites[k] = l_sites[k].lower()
            lige[k] = re.compile(l_sites[k])
        callback = lambda pat: pat.group(0).lower()
        while len(quals) <= nreads:
            try:
                next(fhandler)
            except StopIteration:
                break
            seq = next(fhandler)
            # ligation sites replaced by lower case to ease the search
            for lig in list(l_sites.values()):
                seq = re.sub(lig.upper(), callback, seq)
            for r_enz in r_enzs:
                sites[r_enz].extend([m.start() for m in site[r_enz].finditer(seq)])
                # TODO: you cannot have a repaired/fixed site in the middle of
                # the sequence, this could be only checked at the beginning
                fixes[r_enz].extend([m.start() for m in fixe[r_enz].finditer(seq)])
            for k in lige:  # for each paired of cut-site
                liges[k].extend([m.start() for m in lige[k].finditer(seq)])
                if lige[k].search(seq):
                    ligep[k] += 1
            # store the number of Ns found in the sequences
            if 'N' in seq:
                henes.extend([i for i, s in enumerate(seq) if s == 'N'])
            next(fhandler)
            line = next(fhandler)
            quals.append([phred[i] for i in line.strip()])
    fhandler.close()
    if not nreads:
        nreads = len(quals)
    quals = zip_longest(*quals, fillvalue=float('nan'))
    meanquals, errorquals = list(zip(*[(nanmean(q), nanstd(q)) for q in quals]))
    max_seq_len = len(meanquals)

    if axe:
        ax = axe
        fig = axe.get_figure()
        ax2 = fig.add_subplot(212)
    else:  # configure plot
        if len(r_enzs) == 1 and r_enzs[0] is None:  # do both plots
            _, ax = plt.subplots(1, 1, figsize=(15, 6))
        else:  # only do the quality_plot plot
            _, (ax, ax2) = plt.subplots(2, 1, figsize=(15, 12))
        ax.patch.set_facecolor('lightgrey')
        ax.patch.set_alpha(0.4)
        ax.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
        ax.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
        ax.set_axisbelow(True)
        # remove tick marks
        ax.tick_params(axis='both', direction='out', top=False, right=False,
                       left=False, bottom=False)
        ax.tick_params(axis='both', direction='out', top=False, right=False,
                       left=False, bottom=False, which='minor')

    ax.errorbar(list(range(max_seq_len)), meanquals,
                linewidth=1, elinewidth=1, color='darkblue',
                yerr=errorquals, ecolor='orange')

    ax.set_xlim((0, max_seq_len))
    ax.set_xlabel('Nucleotidic position')
    ax.set_ylabel('PHRED score')
    ax.set_title('Sequencing Quality (%d reads)' % (nreads))
    ax.yaxis.label.set_color('darkblue')
    ax.tick_params(axis='y', colors='darkblue', **tkw)
    axb = ax.twinx()
    # quality_plot plot
    axb.plot([henes.count(i) for i in range(max_seq_len)], linewidth=1,
             color='black', linestyle='--')
    axb.yaxis.label.set_color('black')
    axb.tick_params(axis='y', colors='black', **tkw)
    axb.set_ylabel('Number of "N" per position')
    try: # no Ns found (yes... it happens)
        axb.set_yscale('log')
        with catch_warnings():
            simplefilter("ignore")
            axb.set_ylim((0, axb.get_ylim()[1] * 1000))
    except ValueError:
        axb.set_yscale('linear')
    ax.set_ylim((0, ax.get_ylim()[1]))
    ax.set_xlim((0, max_seq_len))

    # Hi-C plot
    percent_digested_sites = {}
    if not (len(r_enzs) == 1 and r_enzs[0] is None):
        ax.set_title('Sequencing Quality and deconvolution (%s %d reads)' % (
            ', '.join(map(str, r_enzs)), nreads))
        ax.set_xlabel('')
        plt.setp(ax.get_xticklabels(), visible=False)
        ax2.patch.set_facecolor('lightgrey')
        ax2.patch.set_alpha(0.4)
        ax2.grid(ls='-', color='w', lw=1.5, alpha=0.6, which='major')
        ax2.grid(ls='-', color='w', lw=1, alpha=0.3, which='minor')
        ax2.set_axisbelow(True)
        ax2.set_xlabel('Nucleotidic position')

        # seq_len is the length of the line to plot. we don't want to plot
        # if there is no room for the cut-site, or ligation site.
        site_len = max((max([len(r_sites[k]) for k in r_sites]),
                        max([len(l_sites[k]) for k in l_sites]),
                        max([len(d_sites[k]) for k in d_sites])))
        seq_len = max_seq_len - site_len

        # transform dictionaries of positions into dictionaries of counts
        for r_enz in sites:
            sites[r_enz] = [sites[r_enz].count(k) for k in range(seq_len)] # Undigested
            fixes[r_enz] = [fixes[r_enz].count(k) for k in range(seq_len)] # DE
        for r1, r2 in liges:
            liges[(r1, r2)] = [liges[(r1, r2)].count(k) for k in range(seq_len)] # OK

        # in case the pattern of the repaired cut-site contains the target
        # cut-site pattern. These sites were counted twice, once in the
        # undigested, and once in the repaired. We remove them from the
        # repaired:
        for r_enz in r_enzs:
            if d_sites[r_enz] in r_sites[r_enz]:
                pos = r_sites[r_enz].find(d_sites[r_enz])

                fixes[r_enz] = (fixes[r_enz][:pos] +
                                [fixes[r_enz][k] - sites[r_enz][k-pos]
                                 for k in range(pos, seq_len)])
        # same for ligated sites
        for r_enz1 in r_enzs:
            for r_enz2 in r_enzs:
                if d_sites[r_enz1] not in l_sites[(r_enz1, r_enz2)]:
                    continue
                pos = l_sites[(r_enz1, r_enz2)].find(d_sites[r_enz1])
                fixes[r_enz1] = (fixes[r_enz1][:pos] +
                                 [fixes[r_enz1][k] - liges[(r_enz1, r_enz2)][k - pos]
                                  for k in range(pos, seq_len)])

        # remove anything that could be in between the two read ends
        if paired:
            for k in sites:
                sites[k][max_seq_len // 2 - site_len:
                         max_seq_len // 2] = [float('nan')] * site_len
                fixes[k][max_seq_len // 2 - site_len:
                         max_seq_len // 2] = [float('nan')] * site_len
            for k in liges:
                liges[k][max_seq_len // 2 - site_len:
                         max_seq_len // 2] = [float('nan')] * site_len

        # plot undigested cut-sites
        color = iter(plt.cm.Reds(linspace(0.3, 0.95, len(r_enzs))))
        for r_enz in sites:
            # print 'undigested', r_enz
            # print sites[r_enz][:20]
            ax2.plot(sites[r_enz], linewidth=2, color=next(color),
                     alpha=0.9,
                     label='Undigested RE site (%s: %s)' % (r_enz, r_sites[r_enz])
                     if any([f > 0 for f in fixes[r_enz]])
                     else 'Undigested & Dangling-Ends (%s: %s)' % (r_enz, r_sites[r_enz]))
        ax2.set_ylabel('Undigested')
        ax2.yaxis.label.set_color('darkred')
        ax2.tick_params(axis='y', colors='darkred', **tkw)

        lines, labels = ax2.get_legend_handles_labels()

        ax3 = ax2.twinx()
        color = iter(plt.cm.Blues(linspace(0.3, 0.95, len(liges))))
        for r1, r2 in liges:
            # print 'ligated', r1, r2
            # print liges[(r1, r2)][:20]
            ax3.plot(liges[(r1, r2)], linewidth=2, color=next(color),
                     alpha=0.9,
                     label = 'Ligated (%s-%s: %s)' % (r1, r2, l_sites[(r1, r2)].upper()))
        ax3.yaxis.label.set_color('darkblue')
        ax3.tick_params(axis='y', colors='darkblue', **tkw)
        ax3.set_ylabel('Ligated')

        tmp_lines, tmp_labels = ax3.get_legend_handles_labels()
        lines.extend(tmp_lines)
        labels.extend(tmp_labels)

        color = iter(plt.cm.Greens(linspace(0.3, 0.95, len(r_enzs))))
        for i, r_enz in enumerate(r_enzs):
            if any([f > 0 for f in fixes[r_enz]]):
                ax4 = ax2.twinx()
                ax4.spines["right"].set_position(("axes", 1.07))
                make_patch_spines_invisible(ax4)
                ax4.spines["right"].set_visible(True)
                # print 'repaired', r_enz
                # print fixes[r_enz][:20]
                ax4.plot(fixes[r_enz], linewidth=2, color=next(color),
                         alpha=0.9,
                         label='Dangling-ends (%s: %s)' % (r_enz, d_sites[r_enz]))
                ax4.yaxis.label.set_color('darkgreen')
                ax4.tick_params(axis='y', colors='darkgreen', **tkw)
                ax4.set_ylabel('Dangling-ends')
                tmp_lines, tmp_labels = ax4.get_legend_handles_labels()
                lines.extend(tmp_lines)
                labels.extend(tmp_labels)
            else:
                ax2.set_ylabel('Undigested & Dangling-ends')
        ax2.set_xlim((0, max_seq_len))

        # Count ligation sites
        lig_cnt = {}
        for k in liges:
            lig_cnt[k] = (nansum(liges[k]) - liges[k][0] -
                              liges[k][max_seq_len // 2])

        # Count undigested sites
        sit_cnt = {}
        for r_enz in r_enzs:
            sit_cnt[r_enz] = (nansum(sites[r_enz]) - sites[r_enz][0] -
                              sites[r_enz][max_seq_len // 2])

        # Count Dangling-Ends
        des = {}
        for r_enz in r_enzs:
            if any([f > 0 for f in fixes[r_enz]]):
                des[r_enz] = ((100. * (fixes[r_enz][0] + (fixes[r_enz][(max_seq_len // 2)]
                                                          if paired else 0))) / nreads)
            else:
                des[r_enz] = (100. * (sites[r_enz][0] + (sites[r_enz][(max_seq_len // 2)]
                                                         if paired else 0))) / nreads

        # Decorate plot
        title = ''

        for r_enz in r_enzs:
            lcnt = float(sum([lig_cnt[(r_enz1, r_enz2)] * (2 if r_enz1 == r_enz2 else 1)
                              for r_enz1 in r_enzs for r_enz2 in r_enzs
                              if r_enz1 == r_enz or r_enz2 == r_enz]))
            title += ('Percentage of digested sites (not considering Dangling-Ends) '
                      '%s: %.1f%%\n' % (r_enz,
                                        100. * float(lcnt) / (lcnt + sit_cnt[r_enz])))
            percent_digested_sites[r_enz] = 100. * float(lcnt) / (lcnt + sit_cnt[r_enz])
        for r_enz in r_enzs:
            title += 'Percentage of dangling-ends %s: %.1f%%\n' % (r_enz, des[r_enz])

        for r_enz1 in r_enzs:
            for r_enz2 in r_enzs:
                title += ('Percentage of reads with ligation site (%s-%s): %.1f%% \n' %
                          (r_enz1, r_enz2, (ligep[(r_enz1, r_enz2)] * 100.) / nreads))
        plt.title(title.strip(), size=10, ha='left', x=0)
        plt.subplots_adjust(right=0.85)
        ax2.legend(lines, labels, bbox_to_anchor=(0.75, 1.0),
                   loc=3, borderaxespad=0., frameon=False, fontsize=9)
    plt.tight_layout()
    if output_prefix:
        for ext in extension_list:
            plt.savefig(output_prefix + "." + ext, )
        plt.close('all')

    elif not axe:
        plt.show()
    for k in ligep:
        ligep[k] = (ligep[k] * 100.) / nreads
    if len(r_enzs) == 1 and r_enzs[0] is None:
        return {}, {}
    return des, ligep, percent_digested_sites


parser = argparse.ArgumentParser()


parser.add_argument("-f", "--input_fastq", action="store", dest="input_fastq", required=True,
                    help="Input fastq file, can be compressed. Required")
parser.add_argument("-e", "--restriction_enzyme_list", action="store", dest="restriction_enzyme_list", default=None,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of the restriction enzymes.")
parser.add_argument("-n", "--read_number", action="store", dest="read_number", default=100000, type=int,
                    help="Number of reads to analyze. Default: 100000")
parser.add_argument("-s", "--percentage_header_prefix", action="store", dest="percentage_header_prefix", default="",
                    help="Prefix of percentage column in the header. "
                         "Useful if you with to analyze multiple files and combine stats in a single table later. "
                         "Default: '', i.e not set")
parser.add_argument("-i", "--interleaved", action="store_true", dest="interleaved", default=False,
                    help="Input fastq is interleaved, i.e contains both forward and reverse reads. Default: not set")
parser.add_argument("-p", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files.")
parser.add_argument("-x", "--figure_format_list", action="store", dest="figure_format_list", default=["png", "svg"],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of extensions for output figures. Default: png,svg")
args = parser.parse_args()

des, ligep, percent_digested_sites = quality_plot(args.input_fastq, r_enz=args.restriction_enzyme_list,
                                                  nreads=args.read_number,
                                                  axe=None, output_prefix=args.output_prefix,
                                                  extension_list=args.figure_format_list, paired=args.interleaved)

with open(args.output_prefix + ".stats", "w") as out_fd:
    out_fd.write("#type\trestrictase(s)\t{0}percentage\n".format(args.percentage_header_prefix))
    for restrictase in percent_digested_sites:
        out_fd.write("{0}\t{1}\t{2:3.2f}\n".format("digested sites", restrictase, percent_digested_sites[restrictase]))
    for restrictase in des:
        out_fd.write("{0}\t{1}\t{2:3.2f}\n".format("dangling ends", restrictase, des[restrictase]))
    for restrictase_couple in ligep:
        out_fd.write("{0}\t{1}-{2}\t{3:3.2f}\n".format("ligation_sites", restrictase_couple[0],
                                                  restrictase_couple[1], ligep[restrictase_couple]))
