__author__ = 'mahajrod'
import os
from collections import Iterable, OrderedDict
from math import sqrt

import matplotlib.pyplot as plt
import numpy as np

from Collections.GeneralCollections import TwoLvlDict

built_in_flags = {"DA": "desaminase-like",
                  "BR": "location in bad region (masked and so on)",
                  "IP": "indel presence"
                  }


class Record():
    def __init__(self, chrom, pos, description=None, flags=None):
        self.chrom = chrom
        self.pos = pos
        self.description = description
        self.flags = flags

    def __str__(self):
        pass

    def get_location(self, record_dict):
        # function is written for old variant (with sub_feature)s rather then new (with CompoundLocation)
        # id of one SeqRecord in record_dict must be equal to record.pos
        #print(self.description)
        if not self.description:
            self.description = {}

        if "Loc" not in self.description:
            self.description["Loc"] = set([])
        #print(self.chrom, self.pos)
        for feature in record_dict[self.chrom].features:
            if (self.pos - 1) in feature:
                self.description["Loc"].add(feature.type)
                #print(feature)
            for sub_feature in feature.sub_features:
                if (self.pos - 1) in sub_feature:
                    self.description["Loc"].add(sub_feature.type)
        if not self.description["Loc"]:
            self.description["Loc"].add("intergenic")
        #print(self.description)
        #print(self.description["Loc"], self.pos)


class Metadata():
    def __init__(self, metadata=[]):
        self.metadata = metadata

    def add_metadata(self):
        pass

    def __str__(self):
        pass


class Header():
    def __init__(self):
        pass

    def __str__(self):
        pass


class Collection(Iterable):
    # TODO: develop this class to minimize new code in various collections
    def __init__(self, metadata=None, record_list=None, header=None, input_file=None, from_file=False):
        # metadata should be Metadata class
        # header should be Header class
        # record_list should be list of Record class
        if from_file:
            self.read(input_file)
        else:
            self.metadata = metadata
            self.records = record_list
            self.header = header

    def read(self, input_file):
        # collectiontype-dependent function
        pass

    def add_metadata(self):
        # collectiontype-dependent function
        pass

    def add_header(self):
        # collectiontype-dependent function
        pass

    def add_record(self):
        # collectiontype-dependent function
        pass

    def __str__(self):
        collection_string = ""
        if self.metadata:
            collection_string += str(self.metadata)
        if self.header:
            collection_string += "\n" + str(self.header)
        if self.records:
            collection_string += "\n" + "\n".join([str(record) for record in self.records])

        return collection_string

    def __iter__(self):
        for record in self.records:
            yield record

    def filter_records_by_parameter_value(self, parameter, minimum, maximum):
        # TODO: fix
        filtered_records = []
        filtered_out_records = []
        for record in self.records:
            if parameter >= minimum and parameter <= maximum:
                filtered_records.append(record)
            else:
                filtered_out_records.append(record)

    def split_records_by_flags(self, flag_set, mode="all"):
        # possible modes:
        # all - record to be counted as 'with flag' must have all flags from flags_list
        # one - record to be counted as 'with flag' must have at least one flag from flags_list
        flags = set(flag_set)
        with_flags_records = []
        without_flags_records = []
        if mode == "all":
            for record in self:
                if flags & record.flags == flags:
                    with_flags_records.append(record)
                else:
                    without_flags_records.append(record)
        elif mode == "one":
            for record in self:
                if flags & record.flags:
                    with_flags_records.append(record)
                else:
                    without_flags_records.append(record)

        return with_flags_records, without_flags_records

    def write(self, output_file):
        with open(output_file, "w") as out_fd:
            if self.metadata:
                out_fd.write(str(self.metadata) + "\n")
            if self.header:
                out_fd.write(str(self.header) + "\n")
            for record in self.records:
                out_fd.write(str(record) + "\n")

    def get_location(self, record_dict):
        for record in self.records:
            record.get_location(record_dict)

    def _split_regions(self):
        splited_dict = OrderedDict({})
        for record in self.records:
            if record.chrom not in splited_dict:
                splited_dict[record.chrom] = [record]
            else:
                splited_dict[record.chrom].append(record)
        return splited_dict

    def count_locations(self, annotation_black_list=[],
                        allow_several_counts_of_record=False,
                        out_filename="location_counts.t",
                        write=True,
                        count_dir="location_counts"):
        os.system("mkdir -p %s" % count_dir)
        regions_dict = self._split_regions()
        region_counts_dict = TwoLvlDict({})
        for region in regions_dict:
            count_locations_dict = {"intergenic": 0}

            for record in regions_dict[region]:
                if (not record.description["Loc"]) or ("Loc" not in record.description):
                    count_locations_dict["unknown"] += 1
                    continue
                #print(record.description["Loc"])
                if allow_several_counts_of_record:
                    for location in record.description["Loc"]:
                        if location in annotation_black_list:
                            continue
                        if location not in count_locations_dict:
                            count_locations_dict[location] = 1
                        else:
                            count_locations_dict[location] += 1
                else:
                    full_location = []
                    for location in record.description["Loc"]:
                        if location in annotation_black_list:
                            continue
                        full_location.append(location)
                    if not full_location:
                        continue
                    full_location = "/".join(full_location)
                    if full_location not in count_locations_dict:
                        count_locations_dict[full_location] = 1
                    else:
                        count_locations_dict[full_location] += 1

            labels = []
            counts = []
            #colors = []
            for location in count_locations_dict:
                if count_locations_dict[location] == 0 or location in annotation_black_list:
                    continue
                labels.append(location)
                counts.append(count_locations_dict[location])
            region_counts_dict[region] = OrderedDict([(label, count) for label, count in zip(labels, counts)])

        if write:
            region_counts_dict.write("%s/%s" % (count_dir, out_filename))
        return region_counts_dict

    def location_pie(self, pie_name="Location of variants", annotation_colors=[],
                     dpi=150, figsize=(30, 30), facecolor="#D6D6D6",
                     ref_genome=None, explode=True, annotation_black_list=[],
                     allow_several_counts_of_record=False,
                     pie_filename="variant_location_pie.svg",
                     counts_filename="location_counts.t",
                     plot_dir="variant_location_pie",
                     counts_dir="location_counts"):
        print("Drawing location pie...")

        os.system("mkdir -p %s" % plot_dir)
        reference_colors = {"CDS": "#FBFD2B",    # yellow
                            "five_prime_UTR": "#FF000F",     # red
                            "three_prime_UTR": "#000FFF",     # blue
                            "intergenic": "#4ED53F",     # green
                            "other": "#ADB696"
                            }

        if annotation_colors:
            reference_colors = annotation_colors
            reference_colors["other"] = "#ADB696"

        count_locations_dict = self.count_locations(annotation_black_list=annotation_black_list,
                                                    allow_several_counts_of_record=allow_several_counts_of_record,
                                                    out_filename=counts_filename,
                                                    write=True,
                                                    count_dir=counts_dir)
        num_of_regions = len(count_locations_dict)
        side = int(sqrt(num_of_regions))
        if side*side != num_of_regions:
            side += 1
        sub_plot_dict = OrderedDict({})
        fig = plt.figure(2, dpi=dpi, figsize=figsize)
        fig.suptitle(pie_name, fontsize=20)

        index = 1

        for region in count_locations_dict:
            sub_plot_dict[region] = plt.subplot(side, side, index, axisbg="#D6D6D6")

            labels = []
            counts = []
            colors = []
            #print(count_locations_dict)
            for label in count_locations_dict[region]:
                #print(count_locations_dict[region])
                if count_locations_dict[region][label] == 0 or label in annotation_black_list:
                    continue
                labels.append(label)
                counts.append(count_locations_dict[region][label])
                if label not in reference_colors:
                    colors.append(reference_colors["other"])
                else:
                    colors.append(reference_colors[label])

            explodes = np.zeros(len(counts))
            if explode and counts:
                max_count_index = counts.index(max(counts))
                explodes[max_count_index] = 0.1
            plt.pie(counts, explode=explodes, labels=labels, colors=colors,
                    autopct='%1.1f%%', shadow=True, startangle=90)
            plt.title("Region %s" % region, y=1.08)
            # Set aspect ratio to be equal so that pie is drawn as a circle.
            plt.axis('equal')
            index += 1
        plt.savefig("%s/%s" % (plot_dir, pie_filename))
        plt.close()

        """
                print("Drawing location pie...")
        plot_dir = "variant_location_pie"
        os.system("mkdir -p %s" % plot_dir)
        reference_colors = {"CDS": "#FBFD2B",    # yellow
                            "five_prime_UTR": "#FF000F",     # red
                            "three_prime_UTR": "#000FFF",     # blue
                            "intergenic": "#4ED53F",     # green
                            "other": "#ADB696"
                            }

        if annotation_colors:
            reference_colors = annotation_colors
            reference_colors["other"] = "#ADB696"
        count_locations_dict = self.count_locations()
        regions_dict = self._split_regions()
        num_of_regions = len(regions_dict)
        side = int(sqrt(num_of_regions))
        if side*side != num_of_regions:
            side += 1
        sub_plot_dict = OrderedDict({})
        fig = plt.figure(2, dpi=dpi, figsize=figsize)
        fig.suptitle(pie_name, fontsize=20)

        index = 1
        region_counts_dict = OrderedDict({})
        for region in regions_dict:
            sub_plot_dict[region] = plt.subplot(side, side, index, axisbg="#D6D6D6")
            count_locations_dict = {"intergenic": 0}

            for record in regions_dict[region]:
                if (not record.description["Loc"]) or ("Loc" not in record.description):
                    count_locations_dict["intergenic"] += 1
                    continue

                if allow_several_counts_of_record:
                    for location in record.description["Loc"]:
                        if location in annotation_black_list:
                            continue
                        if location not in count_locations_dict:
                            count_locations_dict[location] = 1
                        else:
                            count_locations_dict[location] += 1
                else:
                    full_location = []
                    #print(record.description["Loc"])
                    for location in record.description["Loc"]:
                        if location in annotation_black_list:
                            continue
                        full_location.append(location)
                    full_location = "/".join(full_location)
                    #print(full_location)
                    if full_location not in count_locations_dict:
                        count_locations_dict[full_location] = 1
                    else:
                        count_locations_dict[full_location] += 1
            labels = []
            counts = []
            colors = []
            for location in count_locations_dict:
                if count_locations_dict[location] == 0 or location in annotation_black_list:
                    continue
                labels.append(location)
                counts.append(count_locations_dict[location])
                if location not in reference_colors:
                    colors.append(reference_colors["other"])
                else:
                    colors.append(reference_colors[location])

            explodes = np.zeros(len(counts))
            if explode:
                max_count_index = counts.index(max(counts))
                explodes[max_count_index] = 0.1
            plt.pie(counts, explode=explodes, labels=labels, colors=colors,
                    autopct='%1.1f%%', shadow=True, startangle=90)
            plt.title("Region %s" % region, y=1.08)
            # Set aspect ratio to be equal so that pie is drawn as a circle.
            plt.axis('equal')
            index += 1
            region_counts_dict[region] = dict([(label, count) for label, count in zip(labels, counts)])
        plt.savefig("%s/%s" % (plot_dir, out_file_name))
        plt.close()
        return region_counts_dict
        """

if __name__ == "__main__":
    fg= Collection()