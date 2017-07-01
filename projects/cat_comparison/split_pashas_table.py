__author__ = 'mahajrod'
import os

workdir = "/home/mahajrod/Dropbox/Dobrzhansky/project/cat_comparisons/mobile_elements/"
os.chdir(workdir)
species_list = ["AciJub","CanFam","CarCar","CroCro","FelCat","LynPar","PanLeo","PanOnc","PanPar","PanTig","PriBen","PriViv"]
file_dict = {}
for species in species_list:
    file_dict[species] = open("%s.ltr.bed" % species, "w")

with open("/home/mahajrod/Dropbox/Dobrzhansky/project/cat_comparisons/mobile_elements/AllLtrTable.txt", "r") as in_fd:
    for line in in_fd:
        if line[0] == "#":
            continue
        line_list = line.strip().split("\t")
        species = line_list[0]
        name = line_list[1]
        scaffold, start, end = name.split("_")
        file_dict[species].write("%s\t%s\t%s\t%s\n" % (scaffold, start, end, name))
for species in species_list:
    file_dict[species].close()
