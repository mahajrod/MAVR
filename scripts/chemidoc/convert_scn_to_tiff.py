#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import os
import shutil
import argparse
from RouToolPa.Tools.ImageMagick import Convert
from RouToolPa.Tools.BioFormats import BioConvert
from RouToolPa.Routines import FileRoutines




parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input",
                    type=lambda x: FileRoutines.make_list_of_path_to_files(x.split(",")),
                    help="Comma-separated list of scn files")
parser.add_argument("-t", "--threads", action="store", dest="threads", type=int, default=1,
                    help="Number of threads")

parser.add_argument("-o", "--output_directory", action="store", dest="output", type=FileRoutines.check_path,
                    help="Output directory")
#parser.add_argument("-p", "--convert_options", action="store", dest="convert_options",
#                    help="Options for convert")
parser.add_argument("-d", "--dont_make_negative", action="store_true", dest="dont_negative",
                    help="Dont make negative")

args = parser.parse_args()

temp_dir = "temp/"

FileRoutines.safe_mkdir(temp_dir)

BioConvert.threads = args.threads
Convert.threads = args.threads

BioConvert.parallel_convert(args.input, temp_dir)

if args.dont_negative:
    os.rename(temp_dir, args.output)

else:
    converted_files = os.listdir(temp_dir)
    converted_files = list(map(lambda s: "%s%s" % (temp_dir, s), converted_files))
    Convert.parallel_convert(converted_files, args.output, convert_options=" -negate")#convert_options=args.convert_options)

    shutil.rmtree(temp_dir)
