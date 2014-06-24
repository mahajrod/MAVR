import os


def generate_partition_finder_control_file(alignment_file_name,
                                           genes_file,
                                           output_dir=None,
                                           genes_names_list=None,
                                           branchlengths="linked",
                                           models="raxml",
                                           model_selection="AIC",
                                           search="greedy"):
    #supported models: all, raxml, mrbayes, beast
    #supported model_selection: AIC, AICc, BIC
    #supported search: all, greedy, rcluster, hcluster, user
    print("Generating config file for PartionFinder")
    output_filename = "partition_finder.cfg"
    if output_dir:
        output_filename = output_dir + "/partition_finder.cfg"

    fd_genes = open(genes_file, "r")
    fd_genes.readline()
    coordinates_list = []
    for line in fd_genes:
        line_list = [int(x) for x in line.strip().split("\t")]
        coordinates_list.append(line_list)
    if not line_list[-1]:
        del(line_list[-1])
    fd_genes.close()

    fd = open(output_filename, "w")
    fd.write("""
## ALIGNMENT FILE ##
alignment = %s;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = %s;

## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | beast | <list> ##
##              for PartitionFinderProtein: all_protein | <list> ##
models = %s;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = %s;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]\n
""" % (alignment_file_name, branchlengths, models, model_selection))

    for i in range(0, len(coordinates_list)):
        name = "part%i" % (i + 1)
        if genes_names_list:
            name = genes_names_list[i]
        fd.write("%s_pos1 = %i-%i\\3;\n" % (name, coordinates_list[i][1], coordinates_list[i][2]))
        fd.write("%s_pos2 = %i-%i\\3;\n" % (name, coordinates_list[i][1] + 1, coordinates_list[i][2]))
        fd.write("%s_pos3 = %i-%i\\3;\n" % (name, coordinates_list[i][1] + 2, coordinates_list[i][2]))
    fd.write("""
## SCHEMES, search: all | greedy | rcluster | hcluster | user ##
[schemes]
search = %s;

#user schemes go here if search=user. See manual for how to define.#
""" % search)
    fd.close()


def find_best_partion_model(work_dir, raxml_threads=6):
    os.system("PartitionFinder.py --raxml -p 1 --cmdline-extras '-T %i'  %s" % (raxml_threads, work_dir))


def generate_raxml_partitions_file(partion_finder_output_file, output_filename):
    fd = open(partion_finder_output_file, "r")
    while fd.readline().strip() != "RaxML-style partition definitions":
        pass
    fd_raxml = open(output_filename, "w")
    for line in fd:
        fd_raxml.write(line)
    fd.close()
    fd_raxml.close()


def construct_raxml_tree(alignment_file,
                         partition_file,
                         model="GTRGAMMAI",
                         bootstrap=10000,
                         threads=6,
                         tree_prefix="raxml_tree",
                         path_to_raxml=None):
    raxml_path = "raxml"
    if path_to_raxml:
        raxml_path = path_to_raxml
    os.system("%s -f a -m %s -p 12345 -x 12345 -T %i -s %s -# %i -q %s -n %s" %
              (raxml_path, model, threads, alignment_file, bootstrap, partition_file, tree_prefix))


def add_outgroup(sequences_file, outgroup_file, output_file):
    os.system("cat %s %s > %s" % (sequences_file, outgroup_file, output_file))