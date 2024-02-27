#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from functools import partial
from Data.Nucleotides import degenerated_nucleotides


parser = argparse.ArgumentParser()


parser.add_argument("-r", "--restriction_sites", action="store", dest="restriction_sites", required=True,
                    type=lambda s: s.split(","),
                    help="Comma-separated list of the restriction sites with cutting site labeled. "
                         "For example, '^GATC,G^ANTC'")
parser.add_argument("-s", "--restriction_symbol", action="store", dest="restriction_symbol", default="^",
                    help="Symbol marking the cut position. Default: '^'")
parser.add_argument("-e", "--expand_ambiguous", action="store_true", dest="expand_ambiguous", default=False,
                    help="Expand ambiguous nucleotides. Default: False")

args = parser.parse_args()

left_site_list = []
right_site_list = []
for site in args.restriction_sites:
    if args.restriction_symbol not in site:
        raise ValueError("ERROR! Restriction site '{0}' has no restriction symbol {1}!".format(site,
                                                                                               args.restriction_symbol))
    number_of_cuts = site.count(args.restriction_symbol)
    if number_of_cuts > 1:
        raise ValueError("ERROR! Restriction site '{0}' has more than one cutting site!".format(site))
    site_len = len(site) - 1
    split_list = site.split(args.restriction_symbol)

    left_fragment_len = len(split_list[0])
    right_fragment_len = len(split_list[1])
    full_site = split_list[0] + split_list[1]

    if left_fragment_len > right_fragment_len:
        left_site_list.append(split_list[0])
        right_site_list.append(full_site if right_fragment_len == 0 else full_site[right_fragment_len:])
    elif left_fragment_len == right_fragment_len:
        left_site_list.append(split_list[0])
        right_site_list.append(split_list[1])
    else:
        left_site_list.append(full_site if left_fragment_len == 0 else full_site[:-left_fragment_len])
        right_site_list.append(split_list[1])

ligation_site_list = []
for left_site in left_site_list:
    for right_site in right_site_list:
        ligation_site_list.append(left_site + right_site)


def expand_ambiguous_sites(prefix_list, suffix):
    #print(suffix)
    if suffix == "":
        #print(prefix_list)
        return prefix_list
    new_prefix_list = []
    for prefix in prefix_list:
        for nucleotide in degenerated_nucleotides[suffix[0]]:
            new_prefix_list.append(prefix + nucleotide)
    return expand_ambiguous_sites(new_prefix_list, suffix[1:])


# fast removal of duplicated sites with order remained, works since python 3.7
ligation_site_list = list(dict.fromkeys(ligation_site_list))
#print(ligation_site_list)
if args.expand_ambiguous:
    expanded_ligation_list = []
    for site in ligation_site_list:
        #print(site)
        expanded_ligation_list += expand_ambiguous_sites([""], site)
    #print(expanded_ligation_list)
    expanded_ligation_list = list(dict.fromkeys(expanded_ligation_list))

    for ligation_site in expanded_ligation_list:
        print(ligation_site)
else:
    for ligation_site in ligation_site_list:
        print(ligation_site)

