#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import argparse
from RouToolPa.Routines import ProjectRoutines


parser = argparse.ArgumentParser()

parser.add_argument("-w", "--workdir", action="store", dest="workdir", required=True,
                    help="Working directory")
parser.add_argument("-p", "--project_name", action="store", dest="project_name", required=True,
                    help="Name of the project")
parser.add_argument("-s", "--species_list", action="store", dest="species_list", required=True,
                    type=lambda s: s.split(","),
                    help="SComma-separated list of species")

args = parser.parse_args()

ProjectRoutines.initiate(args.workdir, args.project_name, args.species_list)
