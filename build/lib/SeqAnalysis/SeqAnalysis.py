#!/usr/bin/env python2

import os
import pickle
import matplotlib.pyplot as pyplot
from random import randint

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from Parse.ParseGeneral import parse_alignment


def record_generator(record_dict):
    for record_id in record_dict:
        yield record_dict[record_id]


def record_generator_by_id(record_dict, record_id_list):
    for record_id in record_id_list:
        yield record_dict[record_id]


def get_general_statistics(input_file, file_format, output_file="statistics.t", write_to_file=False):
    record_dict = SeqIO.to_dict(SeqIO.parse(input_file, file_format))
    statistics_dict = {}
    for record_id in record_dict:
        statistics_dict[record_id] = [len(record_dict[record_id].seq), len(record_dict[record_id].features)]
    if write_to_file:
        fd = open(output_file, "w")
        metadata = "#Totaly records\t%i\n" % len(record_dict)
        fd.write(metadata)
        for record_id in sorted(list(record_dict.keys())):
            fd.write("%s\t%i\t%i\n" % (record_id, statistics_dict[record_id][0], statistics_dict[record_id][1]))
        fd.close()
    return statistics_dict


def get_multifile_general_statistics(input_file_list, file_format, output_file="statistics.t", output_file2="file_statistics.t", names_list=None, write_to_file=False):
    statistics_dict = {}
    names_dict = {}
    if names_list:
        for i in range(0, len(input_file_list)):
            names_dict[input_file_list[i]] = names_list[i]
    else:
        for i in range(0, len(input_file_list)):
            names_dict[input_file_list[i]] = input_file_list[i]

    for file_entry in input_file_list:
        print("Gathering statistics from %s" % file_entry)
        statistics_dict[names_dict[file_entry]] = get_general_statistics(file_entry, file_format)

    if write_to_file:
        fd_file = open(output_file2, "w")
        fd = open(output_file, "w")
        metadata = "#Totaly files\t%i\n" % len(input_file_list)
        fd.write(metadata)
        fd_file.write(metadata)
        for entry in sorted(list(statistics_dict.keys())):
            fd.write("%s\t%i\n" % (entry, len(statistics_dict[entry])))
            fd_file.write("%s\t%i\n" % (entry, len(statistics_dict[entry])))
            for record_id in sorted(list(statistics_dict[entry].keys())):
                fd.write("\t%s\t%i\t%i\n" % (record_id, statistics_dict[entry][record_id][0], statistics_dict[entry][record_id][1]))
        fd.close()
        fd_file.close()
    return statistics_dict


def get_statistics(filelist, index_filename, taxonomy_file_prefix="taxonomy", filetype="genbank"):
    record_dict = SeqIO.index_db(index_filename, filelist, filetype)
    number_of_records = len(record_dict)
    print("Gathering statistics...")
    print("Totaly records:\t %i" % number_of_records)
    taxonomy_dict = {}
    #constructing nested taxonomy summary as nested dictionary
    print("Constucting taxonomy distribution...")

    for record_id in record_dict:
        temp_dict = taxonomy_dict
        #print(record_id)
        for taxon in record_dict[record_id].annotations["taxonomy"] + [record_dict[record_id].annotations["organism"]]:
            if taxon not in temp_dict:
                temp_dict[taxon] = {"Total": 1}
            else:
                temp_dict[taxon]["Total"] += 1

            temp_dict = temp_dict[taxon]
        #print(record_dict[record_id].annotations)
    #print(taxonomy_dict)
    with open(taxonomy_file_prefix + ".pickle", 'wb') as fd:
        pickle.dump(taxonomy_dict, fd)
    print_formatted_nested_dict(taxonomy_dict, taxonomy_file_prefix + ".tax", indent_type="\t")

    return record_dict, taxonomy_dict


def get_taxonomy(record_dict, output_prefix="taxonomy"):
    #MODIFIED: to record_dict as input intead of filelist
    #record_dict = SeqIO.index_db(index_filename, filelist, filetype)
    number_of_records = len(record_dict)
    print("Gathering taxonomy statistics...")
    print("Totaly records:\t %i" % number_of_records)
    #species = {}
    genus_dict = {}
    family_dict = {}
    order_dict = {}
    #constructing nested taxonomy summary as nested dictionary
    print("Constucting taxonomy distribution...")

    for record_id in record_dict:
        if record_dict[record_id].annotations["taxonomy"][-1] not in genus_dict:
            genus_dict[record_dict[record_id].annotations["taxonomy"][-1]] = [record_id]
        else:
            genus_dict[record_dict[record_id].annotations["taxonomy"][-1]].append(record_id)

        for taxon in record_dict[record_id].annotations["taxonomy"]:
            if taxon[-6:] == "formes":
                if taxon not in order_dict:
                    order_dict[taxon] = [record_id]
                else:
                    order_dict[taxon].append(record_id)
            if taxon[-4:] == "idae":
                if taxon not in family_dict:
                    family_dict[taxon] = [record_id]
                else:
                    family_dict[taxon].append(record_id)

    def write_taxa_dict(taxa_dict, taxa_file):
        with open(taxa_file, "w") as fd:
            for taxon in taxa_dict:
                fd.write("%s\t%i\t" % (taxon, len(taxa_dict[taxon])) + ",".join(taxa_dict[taxon]) + "\n")

    write_taxa_dict(genus_dict, output_prefix + "_genus.tax")
    write_taxa_dict(family_dict, output_prefix + "_family.tax")
    write_taxa_dict(order_dict, output_prefix + "_order.tax")

    return genus_dict, family_dict, order_dict


def get_taxonomy_from_genbank_files(filelist, index_filename):
    record_dict = SeqIO.index_db(index_filename, filelist, "genbank")
    return get_taxonomy(record_dict)


def count_species(record_dict, output_filename="count_species.count"):
    species_count_dict = {}
    i = 1
    print("Counting species...")
    for record_id in record_dict:
        #print("%i records were handled" % i)
        organism = record_dict[record_id].annotations['organism']
        if organism not in species_count_dict:
            species_count_dict[organism] = [1, [record_id]]
        else:
            species_count_dict[organism][0] += 1
            species_count_dict[organism][1].append(record_id)
        i += 1
    number_of_species = len(species_count_dict)
    print ("Totaly %i species were found" % number_of_species)
    fd = open(output_filename, "w")
    fd.write("#Number of species\t%i\n" % number_of_species)
    for organism in species_count_dict:
        fd.write(organism + "\t%i\t%s\n" % (species_count_dict[organism][0], ",".join(species_count_dict[organism][1])))
        #print(record_dict[record_id].annotations['taxonomy'], record_dict[record_id].annotations['organism'])
    return number_of_species, species_count_dict


def split_records_by_taxa_level(record_dict, prefix, taxa_level=1, filetype="genbank"):
    """taxa levels starts from 0"""
    taxa_dict = {}
    print("Spliting records by %i taxa level" % taxa_level)
    for record_id in record_dict:
        if len(record_dict[record_id].annotations["taxonomy"]) > taxa_level:
            taxon = record_dict[record_id].annotations["taxonomy"][taxa_level]
            if taxon not in taxa_dict:
                taxa_dict[taxon] = [record_id]
            else:
                taxa_dict[taxon].append(record_id)
        else:
            print(record_id,
                  "Not classified for selected taxa level",
                  "Taxonomy:",
                  record_dict[record_id].annotations["taxonomy"],
                  record_dict[record_id].annotations["organism"])

    def taxon_generator():
        for record_id in taxa_dict[taxon]:
            yield record_dict[record_id]

    for taxon in taxa_dict:
        taxon_name = taxon.replace(" ", "_")
        #taxon_records_list = [record_dict[record_id] for record_id in taxa_dict[taxon]]
        #print (taxon_records_list)
        #print(prefix, taxon_name, "%i" % taxa_level)
        filename = str(prefix) + "_" + str(taxa_level) + "_" + str(taxon_name) + ".gb"
        print(filename, len(taxa_dict[taxon]))
        SeqIO.write(taxon_generator(), filename, filetype)

    number_of_taxa = len(taxa_dict)
    print("Totaly %i taxa of %i level were found" % (number_of_taxa, taxa_level))
    return number_of_taxa


def print_formatted_nested_dict(nested_dict, output_file, indent_type="\t"):
    #at moment only \t indent is supported

    def fwrite_nested_dict(nested_dict, fd, indent_type, indent_counter):
        for taxon in nested_dict:
            if taxon != "Total":
                #print(nested_dict[taxon])
                fd.write(indent_type*indent_counter + taxon + "\t%i\n" % nested_dict[taxon]["Total"])
                fwrite_nested_dict(nested_dict[taxon], fd, indent_type, indent_counter + 1)

    fd = open(output_file, "w")
    fwrite_nested_dict(nested_dict, fd, indent_type, 0)
    fd.close()
    return 0


def draw_species_counts_distribution(number_of_species_dict,
                                     distribution_file,
                                     filter_min_value=None,
                                     filter_max_value=None):

    figure = pyplot.figure(dpi=400, figsize=(30, 20))
    axes = figure.add_subplot(1, 1, 1)
    values = []
    print("Drawning distribution of species counts...")
    if filter_min_value and filter_max_value:
        for value in number_of_species_dict.values():
            if value >= filter_min_value  and value <= filter_max_value:
                values.append(value)
        print("Counting only species with sequences in range [%i, %i]" % (filter_min_value, filter_max_value))
    elif filter_min_value:
        for value in number_of_species_dict.values():
            if value >= filter_min_value:
                values.append(value)
        print("Counting only species with % or more sequences" % filter_min_value)
    elif filter_max_value:
        for value in number_of_species_dict.values():
            if value <= filter_max_value:
                values.append(value)
        print("Counting only species with % or less sequences" % filter_max_value)
    else:
        values = list(number_of_species_dict.values())
        print("No filter was set for number of sequences per species")
    max_value = max(values)
    print(values, max_value)
    n, bins, patches = axes.hist(values, max_value, facecolor='green')
    axes.set_xlabel('Number of sequences per species')
    axes.set_ylabel('Number of species')
    axes.set_title('Distribution of sequences per species')
    print("Totaly %i species" % len(values))
    pyplot.savefig(distribution_file)
    axes.set_yscale('log')
    pyplot.savefig('log_' + distribution_file)

    return 0


def filter_by_taxa(record_dict,
                   taxa_list,
                   filter_type="white_list",
                   output_filename="filtered_by_taxa.gb",
                   store_filtered_out=True,
                   filtered_out_filename="filtered_out_by_taxa.gb",
                   output_type="genbank",
                   return_record_generator=False):
    print("Filtering by taxa...")
    print("Totaly %i records" % len(record_dict))


    filtered_id_list = []
    filtered_out_id_list = []
    taxa_set = set(taxa_list)
    if filter_type == "white_list":
        print("Taxa to be retained: %s" % ", ".join(taxa_list))
        for record_id in record_dict:
            #print(record_id)
            if set(record_dict[record_id].annotations["taxonomy"] + [record_dict[record_id].annotations["organism"]]) & taxa_set:
                filtered_id_list.append(record_id)
            else:
                filtered_out_id_list.append(record_id)
    elif filter_type == "black_list":
        print("Taxa to be filtered out: %s" % ", ".join(taxa_list))
        for record_id in record_dict:
            for filter_entry in taxa_set:
                for taxa_entry in record_dict[record_id].annotations["taxonomy"] + [record_dict[record_id].annotations["organism"]]:
                    if filter_entry in taxa_entry:
                        filtered_out_id_list.append(record_id)
                        break
                else:
                    continue
                break
            else:
                filtered_id_list.append(record_id)
    SeqIO.write(record_generator_by_id(record_dict, filtered_id_list), output_filename, output_type)
    if store_filtered_out:
        SeqIO.write(record_generator_by_id(record_dict, filtered_out_id_list), filtered_out_filename, output_type)
    number_of_retained = len(filtered_id_list)
    print("Retained %i records" % number_of_retained)

    if return_record_generator:
        return record_generator_by_id(record_dict, filtered_id_list)


def split_by_taxa(record_dict, taxa_list, output_suffix, output_type="genbank", outfiltered_name="outfiltered"):
    print("Spliting taxa...")
    taxa_dict = {}
    for taxa in taxa_list:
        taxa_dict[taxa] = []
    taxa_dict[outfiltered_name] = []

    for record_id in record_dict:
        record_taxa = record_dict[record_id].annotations["taxonomy"] + [record_dict[record_id].annotations["organism"]]
        for taxa in taxa_list:
            if taxa in record_taxa:
                taxa_dict[taxa].append(record_id)
                break
        else:
            taxa_dict[outfiltered_name].append(record_id)

    if output_type == "genbank":
        out_extension = ".gb"

    for taxa in taxa_dict:
        print("Writing %s taxa" % taxa)
        SeqIO.write(record_generator_by_id(record_dict, taxa_dict[taxa]), taxa + output_suffix + out_extension, output_type)
    return taxa_dict


def filter_by_source(record_dict,
                     source_list,
                     filter_type="white_list",
                     output_filename="filtered_by_source.gb",
                     output_type="genbank",
                     return_record_generator=False):
    print("Filtering by source")
    print("Totaly %i records" % len(record_dict))
    filtered_id_list = []
    filtered_out_id_list = []
    source_set = set(source_list)
    if filter_type == "white_list":
        print("Sources to be retained: %s" % ", ".join(source_list))
        for record_id in record_dict:
            if ("source" not in record_dict[record_id].annotations) or (not (record_dict[record_id].annotations["source"])):
                continue
            #if record_id == "X74571.1" :
            #    continue
            #print(record_id)
            #print(record_dict[record_id])
            """
            products = ""
            for feature in record_dict[record_id].features:
                if "product" in feature.qualifiers:
                    products += " " + feature.qualifiers["product"]
            """
            if set(record_dict[record_id].annotations["source"].split()) & source_set: # or (product_set & source_set):
                filtered_id_list.append(record_id)
            else:
                filtered_out_id_list.append(record_id)
    elif filter_type == "black_list":
        print("Sources to be filtered out: %s" % ", ".join(source_list))
        for record_id in record_dict:
            if "source" not in record_dict[record_id].annotations:
                continue
            """
            product_set = set([])
            for feature in record_dict[record_id].features:
                if "product" in feature.qualifiers:
                    product_set = set(feature.qualifiers["product"])
            """
            if not (set(record_dict[record_id].annotations["source"].split()) & source_set):# or (product_set & source_set)):
                filtered_id_list.append(record_id)
            else:
                filtered_out_id_list.append(record_id)
    SeqIO.write(record_generator_by_id(record_dict, filtered_id_list), output_filename, output_type)
    SeqIO.write(record_generator_by_id(record_dict, filtered_out_id_list), "filtered_out_by_source.gb", output_type)
    number_of_retained = len(filtered_id_list)
    print("Retained %i records" % number_of_retained)

    if return_record_generator:
        return record_generator_by_id(record_dict, filtered_id_list)


def parse_counts_file(species_count_file):
    fd = open(species_count_file, "r")
    number_of_species = fd.readline().strip().split("\t")[-1]
    #print(number_of_species)
    species_count_dict = {}
    for line in fd:
        line_list = line.strip().split("\t")
        #print(line_list)
        species_count_dict[line_list[0]] = [int(line_list[1]), line_list[2].split(",")]
    fd.close()
    return number_of_species, species_count_dict


def sort_by_number(species_count_file, output_filename_prefix, counts_list=[1, 5]):
    number_of_species, species_count_dict = parse_counts_file(species_count_file)
    length_of_counts = len(counts_list)
    fd_list = [open(output_filename_prefix + "_%i_%i_species.count" % (counts_list[i], counts_list[i+1] - 1), "w") for i in range(0, length_of_counts-1)]
    fd_list.append(open(output_filename_prefix + "_%i+_species.count" % (counts_list[-1]), "w"))
    number_of_species = [0 for i in range(0, length_of_counts)]
    species_record_list = [[] for i in range(0, length_of_counts)]

    for species in species_count_dict:
        for i in range(0, length_of_counts-1):
            #print(species_count_dict[species][0])

            if counts_list[i] <= species_count_dict[species][0] < counts_list[i+1]:
                number_of_species[i] += 1
                species_record_list[i].append(species)
                #fd_list[i].write(species + "\t%i\t%s\n" % (species_count_dict[species][0], ",".join(species_count_dict[species][1])))
                break
        else:
            number_of_species[-1] += 1
            species_record_list[-1].append(species)
            #fd_list[-1].write(species + "\t%i\t%s\n" % (species_count_dict[species][0], ",".join(species_count_dict[species][1])))

    for i in range(0, length_of_counts):
        fd_list[i].write("#Number of species\t%i\n" % number_of_species[i])
        for species in species_record_list[i]:
            fd_list[i].write(species + "\t%i\t%s\n" % (species_count_dict[species][0], ",".join(species_count_dict[species][1])))
        fd_list[i].close()

    #for fd in fd_list:
    #    fd.close()


def get_alignment(input_file,
                  output_file="alignment.fasta",
                  aligner_options=None,
                  aligner="MAFFT",
                  path_to_aligner="mafft",
                  threads=6):

    if aligner == "MAFFT":
        if aligner_options:
            aligner_string = "%s --thread %i " % (path_to_aligner, threads) + aligner_options + " " + input_file + " > " + output_file
        else:
            aligner_string = "%s --thread %i --globalpair --maxiterate 1000 " % (path_to_aligner, threads) + input_file + " > " + output_file
            print(aligner_string)
    os.system(aligner_string)


def get_codon_alignment(protein_alignment, nucleotide_seq_dict, codon_alignment_file):
    codon_alignment = {}
    for record in protein_alignment:
        nucleotide_seq = ""
        i = 0
        for aminoacid in record.seq:
            if aminoacid == "-":
                nucleotide_seq += "---"
                continue
            else:
                nucleotide_seq += str(nucleotide_seq_dict[record.id].seq[3*i:3*(i+1)])
                i += 1
        codon_alignment[record.id] = SeqRecord(Seq(nucleotide_seq),
                                               id=record.id,
                                               description=record.description,
                                               name=record.name)
        #print(record.id, record.seq)
    SeqIO.write(list(codon_alignment.values()), codon_alignment_file, "fasta")
    return codon_alignment


def merge_alignment(alignment_file_list, merged_alignment_file):
    print("Merging alignments...")
    alignment_list = []
    sequence_lengthes = []
    #print(alignment_file_list)

    alignment_file_list_sorted = sorted(alignment_file_list)
    print(alignment_file_list_sorted)
    for alignment_file in alignment_file_list_sorted:
        #alignment_file.sort()
        parsed = parse_alignment(alignment_file)
        parsed.sort()
        alignment_list.append(parsed)
    merged_alignment = None
    for alignment in alignment_list:
        if not merged_alignment:
            sequence_lengthes.append(alignment.get_alignment_length())
            merged_alignment = alignment
            continue
        #print(alignment)
        sequence_lengthes.append(alignment.get_alignment_length())
        merged_alignment += alignment
    SeqIO.write(merged_alignment, merged_alignment_file, "fasta")
    sequence_coordinates = []
    #
    for seq_length in sequence_lengthes:
        if not sequence_coordinates:
            sequence_coordinates.append((1, seq_length))
            continue
        sequence_coordinates.append((sequence_coordinates[-1][1]+1, sequence_coordinates[-1][1]+seq_length))
    #print(sequence_coordinates)
    return merged_alignment, sequence_lengthes, sequence_coordinates


def write_sequences_coordinates_to_file(seq_coordinates_tuple_list, output_filename):
    fd = open(output_filename, "w")
    fd.write("#length\tstart\tend\n")
    for coord_tuple in seq_coordinates_tuple_list:
        fd.write("%i\t%i\t%i\n" % (coord_tuple[1] - coord_tuple[0] + 1, coord_tuple[0], coord_tuple[1]))
    fd.close()
    return 1


def test_splited_genes(record_dict, tested_record_dict, output_filename="check_gene_absence.t"):
    #tested_record_dict - key = gene name , element - list of record_ids
    absent_ids_dict = {}
    for key in tested_record_dict:
        absent_ids_dict[key] = set(record_dict.keys()) - set(tested_record_dict[key])
    with open(output_filename, "w") as fd:
        fd.write("#gene\tabsent\tids\n")
        for key in absent_ids_dict:
            fd.write("%s\t%i\t%s\n" % (key, len(absent_ids_dict[key]), ",".join(absent_ids_dict[key])))


def split_mitochondrion_genome_by_genes(record_dict, black_list=[]):
    mithochondrion_synonym_table = [
                                ["12S_rRNA",
                                 "12S rRNA",
                                 "12S ribosomal RNA",
                                 "small subunit ribosomal RNA",
                                 "s-rRNA", "s-rRNA; 12S ribosomal RNA",
                                 "small ribosomal RNA subunit RNA",
                                 "12S ribosomal RNA",
                                 "12S ribosomal RNA subunit",
                                 "12S rivbosomal RNA",
                                 "12S ribosamal RNA",
                                 "l2S ribosomal RNA",
                                 "12 ribosomal RNA",
                                 "12S ribosormal RNA,"
                                 "12 rRNA",
                                 "s-RNA"],
                                ["16S_rRNA",
                                 "16S rRNA",
                                 "16S ribosomal RNA",
                                 "large subunit ribosomal RNA",
                                 "l-rRNA",
                                 "l-rRNA; 16S ribosomal RNA",
                                 "large ribosomal RNA subunit RNA",
                                 "16S ribosomal RNA",
                                 "16S ribosomal RNA subunit",
                                 "16S rivbosomal RNA",
                                 "16S ribosamal RNA",
                                 "l6S ribosomal RNA",
                                 "16 ribosomal RNA",
                                 "16S ribosormal RNA",
                                 "16 rRNA",
                                 "l-RNA"],
                                ["ATP6",
                                 "atp6",
                                 "ATPase6",
                                 "ATPase 6",
                                 "ATPase subunit 6",
                                 "ATP synthase F0 subunit 6",
                                 "ATP synthetase F0 subunit 6",
                                 "ATP synthase subunit 6"
                                 "ATPase subunits 6",
                                 "adenosine triphosphatase subunit 6",
                                 "ATPase subunit-6"],
                                ["ATP8",
                                 "atp8",
                                 "ATPase8",
                                 "ATPase 8",
                                 "ATPase subunit 8",
                                 "ATP synthase F0 subunit 8",
                                 "ATP synthetase F0 subunit 8",
                                 "ATP synthase subunit 8",
                                 "ATPase subunits 8",
                                 "adenosine triphosphatase subunit 8",
                                 "adenosine triphosphate subunit 8",
                                 "ATPase subunit-8"],
                                ["COX1", "COXI",
                                 "cytochrome c oxidase subunit 1",
                                 "cytochrome c oxidase subunit I",
                                 "Cytochrome c oxidase subunit 1",
                                 "cytochrome oxidase subunit I",
                                 "chytochrome c oxidase subunit I",
                                 "COI",
                                 "CO1", "CO 1", "CO I", "coi",
                                 "product: cytochrome c oxidase subunit I",
                                 "cytochrome oxidase subunit 1"],
                                ["COX2", "COXII",
                                 "cytochrome c oxidase subunit 2",
                                 "cytochrome c oxidase subunit II",
                                 "Cytochrome c oxidase subunit 2",
                                 "cytochrome oxidase subunit II",
                                 "chytochrome c oxidase subunit II",
                                 "COII",
                                 "CO2", "CO 2", "CO II", "coii",
                                 "cytochrome oxidase subunit 2"],
                                ["COX3", "COXIII",
                                 "cytochrome c oxidase subunit 3",
                                 "cytochrome c oxidase subunit III",
                                 "Cytochrome c oxidase subunit 3",
                                 "cytochrome oxidase subunit III",
                                 "chytochrome c oxidase subunit III",
                                 "COIII",
                                 "CO3", "CO 3", "CO III", "coiii",
                                 "cytochrome oxidase subunit 3"],
                                ["CYTB", "cytochrome b", "Cytochrome b", "cytb", "Cytb", "Cyt b",
                                 "Cytochrome b apoenzyme",
                                 "cytochrome b apoenzyme",
                                 "cytochrome b; TAA stop codon appears afterpolyadenylation"
                                 ],
                                ["ND1", "nd1", "nd 1", "ND 1", "Nd 1",
                                 "NADH dehydrogenase subunit 1",
                                 "NADH hydrogenase subunit 1",
                                 "subunit 1 of the NADH ubiquinone oxidoreductase complex",
                                 "NADH-1", "NADH1"],
                                ["ND2", "nd2", "nd 2", "ND 2", "Nd 2",
                                 "NADH dehydrogenase subunit 2",
                                 "NADH hydrogenase subunit 2",
                                 "subunit 2 of the NADH ubiquinone oxidoreductase complex",
                                 "#NADH dehydrogenase subunit 2",
                                 "NADH-2", "NADH2"],
                                ["ND3", "nd3", "nd 3", "ND 3", "Nd 3",
                                 "NADH dehydrogenase subunit 3",
                                 "NADH hydrogenase subunit 3",
                                 "subunit 3 of the NADH ubiquinone oxidoreductase complex",
                                 "NADH-3", "NADH3"],
                                ["ND4", "nd4", "nd 4", "ND 4", "Nd 4",
                                 "NADH dehydrogenase subunit 4",
                                 "NADH hydrogenase subunit 4",
                                 "subunit 4 of the NADH ubiquinone oxidoreductase complex",
                                 "NADH-4", "NADH4"],
                                ["ND4L", "nd4l", "nd 4l", "ND 4l", "Nd 4l",
                                 "NADH dehydrogenase subunit 4L",
                                 "NADH hydrogenase subunit 4L",
                                 "NADH-4L", "NADH4L"],
                                ["ND5", "nd5", "nd 5", "ND 5", "Nd 5",
                                 "NADH dehydrogenase subunit 5",
                                 "NADH hydrogenase subunit 5",
                                 "subunit 5 of the NADH ubiquinone oxidoreductase complex",
                                 "NADH-5", "NADH5"],
                                ["ND6", "nd6", "nd 6", "ND 6", "Nd 6",
                                 "NADH dehydrogenase subunit 6",
                                 "NADH hydrogenase subunit 6",
                                 "subunit 6 of the NADH ubiquinone oxidoreductase complex",
                                 "NADH-6", "NADH6",
                                 "NADH dehydrogenase subunit-6"],
                                ["tRNA-Val"], ["tRNA-Leu"], ["tRNA-Phe"], ["tRNA-Pro"], ["tRNA-Thr"],
                                ["tRNA-Glu"], ["tRNA-Ser"], ["tRNA-His"], ["tRNA-Arg"], ["tRNA-Gly"],
                                ["tRNA-Lys"], ["tRNA-Asp"], ["tRNA-Tyr"], ["tRNA-Cys"], ["tRNA-Asn"],
                                ["tRNA-Ala"], ["tRNA-Trp"], ["tRNA-Met"], ["tRNA-Ile"], ["tRNA-Gln"]
                                ]
    protein_gene_list = [
                    "ATP6",
                    "ATP8",
                    "COX1",
                    "COX2",
                    "COX3",
                    "CYTB",
                    "ND1",
                    "ND2",
                    "ND3",
                    "ND4L",
                    "ND4",
                    "ND5",
                    "ND6",
                    ]

    rRNA_gene_list = [
                    "12S_rRNA",
                    "16S_rRNA"
                    ]
    region_dict = {}
    protein_dict = {}

    def check_for_synonym(gene_name, gene_synonym_table):
        for gene_synonym_list in gene_synonym_table:
            if gene_name in gene_synonym_list:
                correct_gene_name = gene_synonym_list[0]
                break
        else:
            correct_gene_name = None
        return correct_gene_name

    def get_feature_gene(feature):
        feature_gene_list = []
        if "product" in feature.qualifiers:
            feature_gene_list.append(check_for_synonym(feature.qualifiers["product"][0], mithochondrion_synonym_table))
        if "gene" in feature.qualifiers:
            feature_gene_list.append(check_for_synonym(feature.qualifiers["gene"][0], mithochondrion_synonym_table))

        if ("gene" not in feature.qualifiers) and ("product" not in feature.qualifiers):
            if "note" in feature.qualifiers:
                feature_gene_list.append(check_for_synonym(feature.qualifiers["note"][0], mithochondrion_synonym_table))
        """
        if "note" in feature.qualifiers:
            feature_gene_list.append(check_for_synonym(feature.qualifiers["note"], mithochondrion_synonym_table))
        """
        for gene_name in feature_gene_list:
            if gene_name:
                return gene_name
        return None

    def modify_record(input_record, record_id, record_name, record_description):
        modified_record = input_record
        modified_record.id = record_id
        modified_record.name = record_name
        modified_record.description = record_description
        return modified_record
    nucleotide_exception_records = []
    protein_exception_records = []
    for record_id in record_dict:
        if record_id in black_list:
            continue
        nucleotide_record_features = {}
        protein_record_features = {}
        for gene in rRNA_gene_list + protein_gene_list:
            nucleotide_record_features[gene] = 0
        for gene in protein_gene_list:
            protein_record_features[gene] = 0
        #print("\t%s" % record_id)
        for feature in record_dict[record_id].features:
            #print("\t%s" % record_id)
            #print(feature)
            if feature.type == "source" or feature.type == "gene" or feature.type == "misc_feature":
                continue
            elif feature.type == "rRNA":
                name = get_feature_gene(feature)
                new_record = modify_record(feature.extract(record_dict[record_id]),
                                           record_id,
                                           name,
                                           "")
                if name:
                    nucleotide_record_features[name] += 1
                    #print(record_id, name, nucleotide_record_features)
                    if name not in region_dict:
                        region_dict[name] = [new_record]
                    else:
                        region_dict[name].append(new_record)
                else:
                    print("Unrecognized rRNA feature\t%s" % (record_id))
            elif feature.type == "tRNA":
                """
                name = get_feature_gene(feature)
                description = ""
                if "note" in feature.qualifiers:
                    description = feature.qualifiers["note"][0]
                new_record = modify_record(feature.extract(record_dict[record_id]),
                                           record_id,
                                           name,
                                           description)
                if name:
                    if name not in region_dict:
                        region_dict[name] = [new_record]
                    else:
                        region_dict[name].append(new_record)
                else:
                    print("Unrecognized tRNA feature\t%s\n %s" % (record_id, str(feature)))
                """
                description = ""
                name = ""
                if "note" in feature.qualifiers:
                    description = feature.qualifiers["note"][0]
                if "product" in feature.qualifiers:
                    name = feature.qualifiers["product"][0]
                new_record = modify_record(feature.extract(record_dict[record_id]),
                                           record_id,
                                           name,
                                           description)
                if "product" in feature.qualifiers:
                    if feature.qualifiers["product"][0] not in region_dict:
                        region_dict[feature.qualifiers["product"][0]] = [new_record]
                    else:
                        region_dict[feature.qualifiers["product"][0]].append(new_record)
                else:
                    print("Unrecognized tRNA feature\t%s" % (record_id))

            elif feature.type == "CDS":
                name = get_feature_gene(feature)
                if name:
                    nucleotide_record_features[name] += 1
                    new_region_record = modify_record(feature.extract(record_dict[record_id]),
                                            record_id,
                                            name,
                                            "")
                    if "translation" not in feature.qualifiers:
                        print("Nontranslated CDS feature\t%s\t%s" % (record_id, str(feature)))
                        continue
                        #print(feature)

                    new_protein_record = SeqRecord(seq=Seq(feature.qualifiers["translation"][0], alphabet=IUPAC.protein),
                                                id=record_id,
                                                name=name,
                                                description="")
                    protein_record_features[name] += 1
                    #print(record_id, feature, new_protein_record)
                    if name not in region_dict:
                        region_dict[name] = [new_region_record]
                        protein_dict[name] = [new_protein_record]
                    else:
                        region_dict[name].append(new_region_record)
                        protein_dict[name].append(new_protein_record)
                else:
                    print("Unrecognized CDS feature\t%s\t%s" % (record_id, str(feature)))
                    #print("\n")
                    #print(feature.qualifiers["product"][0])
                    #print(feature.qualifiers["gene"])
                    #print(feature.qualifiers["product"])
                #print(feature.qualifiers["translation"][0])
                #pass
                #print(feature.qualifiers["gene"])
            elif feature.type == "D-loop":
                new_record = modify_record(feature.extract(record_dict[record_id]),
                                           record_id,
                                           "D-loop",
                                           "")
                if "D-loop" not in region_dict:
                    region_dict["D-loop"] = [new_record]
                else:
                    region_dict["D-loop"].append(new_record)
                #print(feature.extract(record_dict[record_id]))
        for key in nucleotide_record_features:
            if nucleotide_record_features[key] != 1:
                print("Suspicious record %s. Duplication or absense genes" % record_id)
                #print(nucleotide_record_features)
                nucleotide_exception_records.append(record_id)
                break
        for key in protein_record_features:
            if protein_record_features[key] != 1:
                print("Suspicious record %s. Duplication or absense genes" % record_id)
                #print(protein_record_features)
                protein_exception_records.append(record_id)
                break
    nucleotide_exception_records = set(nucleotide_exception_records)
    protein_exception_records = set(protein_exception_records)
    sudpicious_records = nucleotide_exception_records | protein_exception_records
    print("Totally suspicious records : %i" % len(sudpicious_records))
    print("Suspicious records %s: " % ",".join(sudpicious_records))

    with open("suspicious_records.t", "w") as fd_sus:
        fd_sus.write("#id\n" + "\n".join(sudpicious_records))

    def write_data(data_dict, data_type, white_list=None):
        os.system("mkdir -p %s" % data_type)
        os.system("mkdir -p %s/fasta" % data_type)
        os.system("mkdir -p %s/genbank" % data_type)
        for entry in data_dict:
            fd_fasta = open("%s/fasta/" % data_type + entry + ".fasta", "w")
            fd_genbank = open("%s/genbank/" % data_type + entry + ".gb", "w")
            if white_list:
                i = 0
                number_of_records = len(data_dict[entry])
                while i < number_of_records:
                    if data_dict[entry][i].id not in white_list:
                        del(data_dict[entry][i])
                        number_of_records -= 1
                        continue
                    i += 1
            SeqIO.write(data_dict[entry], fd_fasta, "fasta")
            SeqIO.write(data_dict[entry], fd_genbank, "genbank")
            fd_fasta.close()
            fd_genbank.close()

    good_ids = set(record_dict.keys()) - nucleotide_exception_records - protein_exception_records
    """
    regions_set = []
    protein_set = []
    if check_presense_coding_rRNA_genes:
        #regions_list = protein_gene_list + rRNA_gene_list
        for region in protein_gene_list:
            temp_set = []
            for record in region_dict[region]:
                temp_set.append(record.id)
            if not regions_set:
                regions_set = set(temp_set)
                protein_set = set(temp_set)
                continue
            regions_set = regions_set & set(temp_set)
            protein_set = protein_set & set(temp_set)
        for region in rRNA_gene_list:
            for record in region_dict[region]:
                temp_set.append(record.id)
            regions_set = regions_set & set(region_dict[region])
    """
    write_data(region_dict, "nucleotide", white_list=good_ids)
    write_data(protein_dict, "protein", white_list=good_ids)
    tested_nucleotide_dict = {}
    tested_protein_dict = {}
    for gene in protein_gene_list + rRNA_gene_list:
        try:
            tested_nucleotide_dict[gene] = list(SeqIO.index("nucleotide/genbank/%s.gb" % gene, "genbank").keys())
        except ValueError:
            print("Error. Possible duplication of %s" % gene)
    for gene in protein_gene_list:
        try:
            tested_protein_dict[gene] = list(SeqIO.index("protein/genbank/%s.gb" % gene, "genbank").keys())
        except ValueError:
            print("Error. Possible duplication of %s in one of genomes" % gene)

    test_splited_genes(record_dict, tested_nucleotide_dict, output_filename="nucleotide_check_absence.t")
    test_splited_genes(record_dict, tested_protein_dict, output_filename="protein_check_absence.t")

"""
def parse_count_species_file(count_species_file):
    fd = open(count_species_file, "r")
    number_of_species = fd.readline().strip().split("\t")[1]
    print(number_of_species)
    for line in fd:
    fd.close()
"""


def get_random_species_genomes(record_dict,
                               count_species_file,
                               output_file,
                               selected_species_file="selected_species.t",
                               output_type="fasta",
                               prev_id_dict={}):
    print("Extracting random genomes(one per species)...")
    #return ids of random genomes, one per species
    fd = open(count_species_file, "r")
    fd_species = open(selected_species_file, "w")
    fd_species.write("#species\tid\n")
    number_of_species = int(fd.readline().strip().split("\t")[1])
    print("Totaly %s species " % number_of_species)
    random_ids = []
    retained_ids = []
    for line in fd:
        line_list = line.strip().split("\t")
        species = line_list[0]
        if species in prev_id_dict:
            if prev_id_dict[species] in record_dict:
                new_id = prev_id_dict[species]
                retained_ids.append(prev_id_dict[species])
        else:
            num_of_genomes = int(line_list[1])
            id_list = line_list[2].split(",")
            random_number = randint(1, num_of_genomes)
            new_id = id_list[random_number-1]
        random_ids.append(new_id)
        fd_species.write("%s\t%s\n" % (species, new_id))
    fd.close()
    fd_species.close()
    random_ids = set(random_ids)
    SeqIO.write([record_dict[id_entry] for id_entry in random_ids], open(output_file, "w"), output_type)
    #print(len(retained_ids), retained_ids)
    return random_ids


def slice_alignment(input_file, output_file, start, end):
    alignment = parse_alignment(input_file)
    SeqIO.write(alignment[start:end], output_file, "fasta")