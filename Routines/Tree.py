__author__ = 'mahajrod'

from ete2 import Tree


class TreeRoutines:

    @staticmethod
    def unroot_tree_from_file(input_tree_file, output_tree_file, input_tree_format, output_tree_format=None):
        with open(input_tree_file, "r") as in_fd:
            with open(output_tree_file, "w") as out_fd:
                for line in in_fd:
                    tree_line = line.strip()
                    tree = Tree(tree_line, format=input_tree_format)
                    tree.unroot()
                    out_fd.write(tree.write(format=output_tree_format if output_tree_format
                                                                    else input_tree_format) + "\n")
