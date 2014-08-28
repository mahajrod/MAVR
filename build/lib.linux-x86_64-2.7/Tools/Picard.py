__author__ = 'mahajrod'
import os
from Routines.Functions import check_path
from Tools.Abstract import Tool

def make_fasta_dict(fasta_file, dict_name, PICARD_dir=""):
    picard_dir = check_path(PICARD_dir)
    os.system("java -jar %sCreateSequenceDictionary.jar R= %s O= %s"
              % (picard_dir, fasta_file, dict_name))


def add_header2bam(input_bam, output_bam, RGID, RGLB, RGPL, RGSM, RGPU, PICARD_dir=""):
    picard_dir = check_path(PICARD_dir)
    os.system("java -XX:MaxDirectMemorySize=4G -jar %sAddOrReplaceReadGroups.jar I= %s O= %s SORT_ORDER=coordinate RGID=%s RGLB=%s  RGPL=%s RGSM=%s RGPU=%s CREATE_INDEX=True"
              % (picard_dir, input_bam, output_bam, RGID, RGLB, RGPL, RGSM, RGPU))

if __name__ == "__main__":
    workdir = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.6m"
    fasta_file = "LAN210_v0.6m.fasta"
    dict_name = "LAN210_v0.6m.dict"
    os.chdir(workdir)
    make_fasta_dict(fasta_file, dict_name,
                    PICARD_dir="/home/mahajrod/Repositories/genetic/NGS_tools/picard-tools-1.115/picard-tools-1.115")
