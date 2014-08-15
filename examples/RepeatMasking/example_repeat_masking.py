#!/usr/bin/env python
import os
from Tools.RepeatSearchTools import RepeatModeler_search, RepeatMasker_search, TRF_search, extract_repbase


#reference_name = "LAN210_v0.7m"
#reference_name = "Alsu24m_no_ambigious"
reference_name = "LAN210_v0.9m"
reference_file = reference_name + ".fasta"
#workdir = "/home/mahajrod/genetics/desaminases/data/%s/masking" % reference_name
#workdir = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24m/masking"
workdir = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.9m/masking"
os.system("mkdir -p %s" % workdir)
os.chdir(workdir)

os.system("ln -fs ../%s %s" % (reference_file, reference_file))
os.system("mkdir -p repeatmodeler custom_lib")
custom_lib_dir = workdir + "/custom_lib"

os.chdir("repeatmodeler")

RepeatModeler_search("../%s" % reference_file, reference_name,
                     RepeatModeler_dir="/home/mahajrod/Repositories/genetic/NGS_tools/RepeatModeler")
os.chdir(workdir)

extract_repbase("fungi", output_file="custom_lib/RepBase_fungi.fasta",
                RepeatMaskerUtils_dir="/home/mahajrod/Repositories/genetic/NGS_tools/RepeatMasker/util")


for entry in os.listdir(workdir + "/repeatmodeler"):
    if "RM" in entry:
        repeatmodeler_dir = workdir + "/repeatmodeler/" + entry
        break

os.chdir(workdir)
os.system("/bin/cp -f %s/consensi.fa.classified %s/custom_lib.fasta" % (repeatmodeler_dir, custom_lib_dir))
os.system("tail -n +2 %s/RepBase_fungi.fasta >> %s/custom_lib.fasta" % (custom_lib_dir, custom_lib_dir))

os.chdir(workdir)
os.system("mkdir -p repeatmasker")
os.chdir("repeatmasker")
os.system("ln -fs ../%s %s" % (reference_file, reference_file))
os.system("ln -fs ../custom_lib/custom_lib.fasta custom_lib.fasta")
RepeatMasker_search(reference_file, "fungi", custom_lib_path="custom_lib.fasta")

os.chdir(workdir)
os.system("mkdir -p TRF")
os.chdir("TRF")
os.system("ln -fs ../repeatmasker/%s.masked %s_masked_repeatmasker.fasta" % (reference_file, reference_name))
TRF_search("%s_masked_repeatmasker.fasta" % reference_name)
os.system("/bin/cp -f *.mask ../%s_masked.fasta" % reference_name)