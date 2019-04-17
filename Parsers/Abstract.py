__author__ = 'mahajrod'
import os
from collections import Iterable, OrderedDict
from math import sqrt

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.ioff()
import numpy as np

from CustomCollections.GeneralCollections import TwoLvlDict

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

    @staticmethod
    def get_synonym(name, use_synonym=False, synonym_dict=None):
        if (not use_synonym) or (not synonym_dict):
            return name

        if name not in synonym_dict:
            return name
        return synonym_dict[name]

    def get_location(self, record_dict, key="Loc", strand_key="strand", feature_type_black_list=[],
                    use_synonym=False, synonym_dict=None):
        # function is written for old variant (with sub_feature)s rather then new (with CompoundLocation)
        # id of one SeqRecord in record_dict must be equal to record.pos
        # locations will be written to description dictionary of record using "key" as key 
        if not self.description:
            self.description = {}

        if key not in self.description:
            self.description[key] = set([])

        if strand_key not in self.description:
            # by default
            self.description[strand_key] = None
        #print(self.chrom, self.pos)
        for feature in record_dict[self.chrom].features:
            if (self.pos - 1) in feature:
                self.description[key].add(self.get_synonym(feature.type, use_synonym=use_synonym,
                                                           synonym_dict=synonym_dict))
               MACE if self.description[strand_key] is None:
                    self.description[strand_key] = feature.strand
                elif feature.strand != self.description[strand_key]:
                    self.description[strand_key] = 0

                #print(feature)
            for sub_feature in feature.sub_features:
                if (self.pos - 1) in sub_feature:
                    self.description[key].add(self.get_synonym(sub_feature.type, use_synonym=use_synonym,
                                                               synonym_dict=synonym_dict))
                    if self.description[strand_key] is None:
                        self.description[strand_key] = sub_feature.strand
                    elif sub_feature.strand != self.description[strand_key]:
                        self.description[strand_key] = 0

        if not self.description[key]:
            # igc == intergenic
            self.description[key].add("igc")

    def check_location(self, bad_region_collection_gff, expression="bad_region.start <= self.pos <= bad_region.end"):
        if not self.flags:
            self.flags = set()
        for bad_region in bad_region_collection_gff:
            if self.chrom != bad_region.chrom:
                continue
            if eval(expression):
                self.flags.add("BR")


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

    def __getitem__(self, item):
        #string-like access
        return self.records[item]

    def __len__(self):
        return len(self.records)

    def pop(self, index=None):
        # i am a fairy idiot - forgot that 0 is equal to False in python.
        # What to do? Everywhere make direct check if variable is None
        return self.records.pop() if index is None else self.records.pop(index)

    def filter_records_by_parameter_value(self, parameter, minimum=None, maximum=None):
        # TODO: check
        filtered_records = []
        filtered_out_records = []
        if (minimum is not None) and (maximum is not None):
            comparison = lambda par, min_value, max_value: min_value <= par <= max_value
        elif (minimum is None) and (maximum is None):
            raise ValueError("Both minimum and maximum thresholds were not set")
        elif minimum is None:
            comparison = lambda par, min_value, max_value: par <= max_value
        else:
            comparison = lambda par, min_value, max_value: min_value <= par

        for record in self.records:
            if comparison(getattr(record, parameter), minimum, maximum):
                filtered_records.append(record)
            else:
                filtered_out_records.append(record)

        return filtered_records, filtered_out_records

    def filter_records_by_expression(self, expression):
        filtered_records = []
        filtered_out_records = []
        for record in self.records:
            if eval(expression):
                filtered_records.append(record)
            else:
                filtered_out_records.append(record)

        return filtered_records, filtered_out_records

    def filter_records_by_function(self, function):
        filtered_records = []
        filtered_out_records = []
        for record in self.records:
            #print("a\na\na\na\na\na\na\na")
            #print(record.description["Power"])
            if function(record):
                filtered_records.append(record)
            else:
                filtered_out_records.append(record)

        return filtered_records, filtered_out_records

    def check_location(self, bad_region_collection_gff, expression="bad_region.start <= self.pos <= bad_region.end"):
        for record in self:
            record.check_location(bad_region_collection_gff, expression=expression)
    """
    def filter_by_expression(self, expression):
        self_type = self.__class__.__name__
        filtered_records, filtered_out_records = self.filter_records_by_expression(expression)
        #TODO: think how to rewrite following damned string!!!!!!!!!!
        exec("from Parsers.%s import %s" % (self_type[10:], self_type))
        return eval(self_type)(metadata=self.metadata, record_list=filtered_records,
                               header_list=self.header_list, from_file=False), \
               eval(self_type)(metadata=self.metadata, record_list=filtered_out_records,
                               header_list=self.header_list, from_file=False)
    """
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

    def get_location(self, record_dict, key="Loc", use_synonym=False, synonym_dict=None):
        for record in self.records:
            record.get_location(record_dict, key=key, use_synonym=use_synonym, synonym_dict=synonym_dict)

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
            count_locations_dict = {"igc": 0}

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
                    full_location.sort()
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
                     full_genome_pie_filename="variant_location_pie_full_genome.svg",
                     counts_filename="location_counts.t",
                     plot_dir="variant_location_pie",
                     counts_dir="location_counts",
                     radius = 1):
        print("Drawing location pie...")

        os.system("mkdir -p %s" % plot_dir)
        reference_colors = {"CDS": "#FBFD2B",    # yellow
                            "5'_UTR": "#FF000F",
                            "five_prime_UTR": "#FF000F",     # red
                            "3'_UTR": "#000FFF",
                            "three_prime_UTR": "#000FFF",     # blue
                            "igc": "#4ED53F",     # green
                            "ncRNA": 'cyan',
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
        all_labels = []
        all_counts = []
        all_colors = []

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
            for i in range(0, len(labels)):
                if labels[i] not in all_labels:
                    all_labels.append(labels[i])
                    all_counts.append(counts[i])
                    all_colors.append(colors[i])
                else:
                    label_index = all_labels.index(labels[i])
                    all_counts[label_index] += counts[i]

        plt.savefig("%s/%s" % (plot_dir, pie_filename))
        plt.close()

        fig = plt.figure(3, dpi=dpi, figsize=(6, 6))
        fig.suptitle(pie_name, fontsize=20)
        plt.subplot(1, 1, 1, axisbg="#D6D6D6")
        all_explodes = np.zeros(len(all_counts))
        if explode and all_counts:
            max_count_index = all_counts.index(max(all_counts))
            all_explodes[max_count_index] = 0.1
        if len(all_labels) > 0:
            max_label_length = max([len(x) for x in all_labels])
            max_count = max(all_counts)
            max_letters = 1
            while int(max_count / 10**max_letters) != 0:
                max_letters += 1

        patches, texts = plt.pie(all_counts, explode=all_explodes, colors=all_colors,
                                            shadow=True, startangle=90, radius=2) #labels=all_labels,

        #patches, texts, autotexts = plt.pie(all_counts, explode=all_explodes, colors=all_colors,
        #                                    shadow=True, startangle=90, autopct='%1.1f%%', radius=4) #labels=all_labels

        #colors = ['yellowgreen','red','gold','lightskyblue','white','lightcoral','blue','pink', 'darkgreen','yellow','grey','violet','magenta','cyan']
        porcent = 100 * np.array(all_counts).astype(np.float32, copy=False)/sum(all_counts)

        labels = ['{0}  -  {1} ({2:1.2f}%)'.format(i.ljust(max_label_length), str(j).ljust(max_letters), k) for i, j, k in zip(all_labels, all_counts, porcent)]

        if len(all_labels) > 0:
            sort_legend = True
            if sort_legend:
                patches, labels, dummy = zip(*sorted(zip(patches, labels, all_counts),
                                                     key=lambda x: x[2],
                                                     reverse=True))

        plt.legend(patches, labels, loc='center left',  fontsize=13, bbox_to_anchor=(-1.0, 0.5)) # bbox_to_anchor=(-0.1, 1.),
        plt.title("Full genome")
            # Set aspect ratio to be equal so that pie is drawn as a circle.
        plt.axis('equal')
        plt.savefig("%s/%s" % (plot_dir, full_genome_pie_filename), bbox_inches='tight')
        plt.close()

if __name__ == "__main__":
    fg = Collection()