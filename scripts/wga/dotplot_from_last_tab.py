#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import argparse
from RouToolPa.Parsers.LAST import CollectionLast
from RouToolPa.Collections.General import IdList, SynDict
from RouToolPa.Routines import DrawingRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input_last_tab", action="store", dest="input_last_tab", required=True,
                    help="File with LAST output in tab format")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")
parser.add_argument("--format", action="store", dest="format", default="tab",
                    help="Format of LAST alignments. Allowed: tab(default), tab_mismap")
parser.add_argument("-w", "--white_target_id_file", action="store", dest="white_target_id_file",
                    help="File with target scaffold ids from white list or corresponding comma-separated list."
                         "NOTE: filtering is done BEFORE renaming by synonyms!")
parser.add_argument("-b", "--black_target_id_file", action="store", dest="black_target_id_file",
                    help="File with target scaffold ids from black list or corresponding comma-separated list."
                         "NOTE: filtering is done BEFORE renaming by synonyms!")

parser.add_argument("-x", "--white_query_id_file", action="store", dest="white_query_id_file",
                    help="File with query scaffold ids from white list or corresponding comma-separated list."
                         "NOTE: filtering is done BEFORE renaming by synonyms!")
parser.add_argument("-c", "--black_query_id_file", action="store", dest="black_query_id_file",
                    help="File with query scaffold ids from black list or corresponding comma-separated list."
                         "NOTE: filtering is done BEFORE renaming by synonyms!")

parser.add_argument("-s", "--query_syn_file", action="store", dest="query_syn_file",
                    help="File with query scaffold id synonyms")
parser.add_argument("--query_syn_file_key_column", action="store", dest="query_syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in query synonym file")
parser.add_argument("--query_syn_file_value_column", action="store", dest="query_syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in query synonym file synonym")

parser.add_argument("-y", "--target_syn_file", action="store", dest="target_syn_file",
                    help="File with target scaffold id synonyms")
parser.add_argument("--target_syn_file_key_column", action="store", dest="target_syn_file_key_column",
                    default=0, type=int,
                    help="Column(0-based) with key(current id) for scaffolds in target synonym file")
parser.add_argument("--target_syn_file_value_column", action="store", dest="target_syn_file_value_column",
                    default=1, type=int,
                    help="Column(0-based) with value(synonym id) for scaffolds in target synonym file synonym")

parser.add_argument("-u", "--target_order_file", action="store", dest="target_order_file",
                    help="File with order of target scaffolds on the plot or corresponding comma-separated list."
                         "NOTE: ordering is done AFTER renaming by synonyms!")
parser.add_argument("-z", "--query_order_file", action="store", dest="query_order_file",
                    help="File with order of query scaffolds on the plot or corresponding comma-separated list."
                         "NOTE: ordering is done AFTER renaming by synonyms!")

parser.add_argument("-l", "--target_label", action="store", dest="target_label",
                    help="Label for target genome(X axis)")
parser.add_argument("-r", "--query_label", action="store", dest="query_label",
                    help="Label for query genome(Y axis)")
parser.add_argument("-t", "--title", action="store", dest="title",
                    help="Title of dot plot")
parser.add_argument("-e", "--extensions", action="store", dest="extensions", type=lambda x: x.split(","),
                    default=["png", ],
                    help="Comma-separated list of extensions for histogram files. Default: png")

parser.add_argument("-d", "--dpi", action="store", dest="dpi", type=int, default=400,
                    help="DPI of figure. Default: 400")
parser.add_argument("-f", "--figsize", action="store", dest="figsize", type=lambda s: map(int, s.split(",")),
                    default=(12, 12),
                    help="Size of figure in inches(two comma-separated ints). Default: 12,12")
parser.add_argument("-a", "--antialiasing", action="store_true", dest="antialiasing", default=False,
                    help="Enable antialiasing. Use this option only for small sequences, i.e segments of chromosomes ")

parser.add_argument("--linewidth", action="store", dest="linewidth", type=float,
                    default=0.01,
                    help="Width of alignment lines. Default: 0.01")
parser.add_argument("--gridwidth", action="store", dest="gridwidth", type=float,
                    default=1,
                    help="Width of grid lines. Default: 1")
parser.add_argument("--scaffold_label_fontsize", action="store", dest="scaffold_label_fontsize", type=float,
                    default=13,
                    help="Fontsize for scaffold labels. Default: 13")
parser.add_argument("--grid_color", action="store", dest="grid_color", default='black',
                    help="Color of grid lines. Default: 'black'")
parser.add_argument("--hide_grid", action="store_true", dest="hide_grid", default=False,
                    help="Hide grid. Default: False")
parser.add_argument("--bar_color", action="store", dest="bar_color", default='grey',
                    help="Color of bars. Default: 'grey'")
parser.add_argument("--same_strand_color", action="store", dest="same_strand_color", default='blue',
                    help="Color of alignment line if query and target are in same strand. Default: 'blue'")
parser.add_argument("--diff_strand_color", action="store", dest="diff_strand_color", default='red',
                    help="Color of alignment line if query and target are in different strand. Default: 'red'")

parser.add_argument("--target_scaffold_labels_angle", action="store", dest="target_scaffold_labels_angle",
                    default=45, type=int,
                    help="Angle for labels of target scaffolds. Default: 45")
parser.add_argument("--query_scaffold_labels_angle", action="store", dest="query_scaffold_labels_angle",
                    default=0, type=int,
                    help="Angle for labels of query scaffolds. Default: 0")
parser.add_argument("--hide_target_labels", action="store_true", dest="hide_target_labels", default=False,
                    help="Hide labels of target scaffolds. Default: False")
parser.add_argument("--hide_query_labels", action="store_true", dest="hide_query_labels", default=False,
                    help="Hide labels of query scaffolds. Default: False")
parser.add_argument("--bottom_offset", action="store", dest="bottom_offset",
                    default=0.1, type=float,
                    help="Bottom offset for subplot. Default: 0.1")
parser.add_argument("--top_offset", action="store", dest="top_offset",
                    default=0.9, type=float,
                    help="Top offset for subplot. Default: 0.9")
parser.add_argument("--left_offset", action="store", dest="left_offset",
                    default=0.1, type=float,
                    help="Left offset for subplot. Default: 0.1")
parser.add_argument("--right_offset", action="store", dest="right_offset",
                    default=0.9, type=float,
                    help="Right offset for subplot. Default: 0.9")
parser.add_argument("--x_axis_visible", action="store_true", dest="x_axis_visible",
                    default=False,
                    help="Make X axis visible. Default: False")
parser.add_argument("--y_axis_visible", action="store_true", dest="y_axis_visible",
                    default=False,
                    help="Make Y axis visible. Default: False")

args = parser.parse_args()

if args.white_target_id_file:
    target_white_list = IdList(filename=args.white_target_id_file) if os.path.isfile(args.white_target_id_file) else IdList(args.white_target_id_file.split(","))
else:
    target_white_list = IdList()

if args.black_target_id_file:
    target_black_list = IdList(filename=args.black_target_id_file) if os.path.isfile(args.black_target_id_file) else IdList(args.black_target_id_file.split(","))
else:
    target_black_list = IdList()

if args.white_query_id_file:
    query_white_list = IdList(filename=args.white_query_id_file) if os.path.isfile(args.white_query_id_file) else IdList(args.white_query_id_file.split(","))
    #print query_white_list
else:
    query_white_list = IdList()

if args.black_query_id_file:
    query_black_list = IdList(filename=args.black_query_id_file) if os.path.isfile(args.black_query_id_file) else IdList(args.black_query_id_file.split(","))
else:
    query_black_list = IdList()

query_syn_dict = SynDict(filename=args.query_syn_file,
                         key_index=args.query_syn_file_key_column,
                         value_index=args.query_syn_file_value_column)
target_syn_dict = SynDict(filename=args.target_syn_file,
                          key_index=args.target_syn_file_key_column,
                          value_index=args.target_syn_file_value_column)

if args.target_order_file:
    target_order_list = IdList(filename=args.target_order_file) if os.path.isfile(args.target_order_file) else IdList(args.target_order_file.split(","))
else:
    target_order_list = IdList()

if args.query_order_file:
    query_order_list = IdList(filename=args.query_order_file) if os.path.isfile(args.query_order_file) else IdList(args.query_order_file.split(","))
else:
    query_order_list = IdList()

last_collection = CollectionLast(args.input_last_tab,
                                 target_white_list=target_white_list,
                                 target_black_list=target_black_list,
                                 query_white_list=query_white_list,
                                 query_black_list=query_black_list,
                                 query_syn_dict=query_syn_dict,
                                 target_syn_dict=target_syn_dict,
                                 format=args.format
                                 )

last_collection.write("%s.syn.tab" % args.output_prefix)
print args.figsize
DrawingRoutines.draw_dot_plot_from_last_alignment(last_collection,
                                                  output_prefix=args.output_prefix,
                                                  extension_list=args.extensions,
                                                  target_black_list=(),
                                                  target_white_list=(),
                                                  target_ordered_list=target_order_list,
                                                  #target_reverse_list=(),
                                                  query_black_list=(),
                                                  query_white_list=(),
                                                  query_ordered_list=query_order_list,
                                                  #query_reverse_list=(),
                                                  figsize=args.figsize, dpi=args.dpi,
                                                  grid_color=args.grid_color,
                                                  bar_color=args.bar_color,
                                                  same_strand_color=args.same_strand_color,
                                                  diff_strand_color=args.diff_strand_color,
                                                  title=args.title,
                                                  target_label=args.target_label,
                                                  query_label=args.query_label,
                                                  linewidth=args.linewidth,
                                                  gridwidth=args.gridwidth,
                                                  antialiased_lines=args.antialiasing,
                                                  scaffold_label_fontsize=args.scaffold_label_fontsize,
                                                  target_scaffold_labels_angle=args.target_scaffold_labels_angle,
                                                  query_scaffold_labels_angle=args.query_scaffold_labels_angle,
                                                  show_grid=not args.hide_grid,
                                                  show_target_labels=not args.hide_target_labels,
                                                  show_query_labels=not args.hide_query_labels,
                                                  top_offset=args.top_offset,
                                                  bottom_offset=args.bottom_offset,
                                                  left_offset=args.left_offset,
                                                  right_offset=args.right_offset,
                                                  x_axis_visible=args.x_axis_visible,
                                                  y_axis_visible=args.y_axis_visible
                                                  )
