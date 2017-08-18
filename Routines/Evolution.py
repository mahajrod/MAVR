__author__ = 'mahajrod'

import os
from collections import OrderedDict
from CustomCollections.GeneralCollections import IdList
from Routines.Phylogenetics import PhylogeneticsRoutines


class EvolutionRoutines(PhylogeneticsRoutines):
    def __init__(self):
        PhylogeneticsRoutines.__init__(self)

    def split_paml_bootstrap_samples(self, bootstrap_file, output_directory, output_prefix,
                                     output_extension="phy"):
        self.safe_mkdir(output_directory)
        counter = 1
        with open(bootstrap_file, "r") as bootstrap_fd:
            for line in bootstrap_fd:
                #print line
                seq_number, pos_number = map(int, line.strip().split())

                output_file = "%s/%s_%i.%s" % (output_directory,
                                               output_prefix,
                                               counter,
                                               output_extension)
                with open(output_file, "w") as out_fd:
                    out_fd.write(line)
                    for i in range(0, seq_number):
                        out_fd.write(bootstrap_fd.next())
                bootstrap_fd.next()  # read empty line between samples
                counter += 1

    def combine_ds_dn_w_from_bootstrap_data(self, input_dir, output_dir, use_node_names_if_possible=True):

        dn_dir = "%s/dN/" % output_dir
        ds_dir = "%s/dS/" % output_dir
        w_dir = "%s/W/" % output_dir

        for directory in output_dir, dn_dir, ds_dir, w_dir:
            self.safe_mkdir(directory)

        input_files = map(lambda s: "%s/%s" % (input_dir, s), os.listdir(input_dir))

        data_dict = OrderedDict()
        for filename in input_files:
            with open(filename, "r") as in_fd:
                in_fd.readline()  # read header
                for line in in_fd:
                    node_id, node_name, dn, ds, w = line.strip().split("\t")

                    if use_node_names_if_possible:
                        node = node_id if node_name == "." else node_name
                    else:
                        node = node_id

                    if node not in data_dict:
                        data_dict[node] = OrderedDict()
                        for parameter in "dN", "dS", "W":
                            data_dict[node][parameter] = IdList()
                    data_dict[node]["dN"].append(dn)
                    data_dict[node]["dS"].append(ds)
                    data_dict[node]["W"].append(w)

        for node in data_dict:
            for parameter in "dN", "dS", "W":
                out_file = "%s/%s/%s.tsv" % (output_dir, parameter, node)
                data_dict[node][parameter].write(out_file)
