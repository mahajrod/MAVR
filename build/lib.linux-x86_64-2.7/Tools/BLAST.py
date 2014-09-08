__author__ = 'mahajrod'
import os
from Tools.Abstract import Tool

class BLAST_plus():

    def make_blast_plus_db(self, input_file, mask_output_file, db_name):
        #makes BLAST database from fasta file
        os.system("dustmasker -in %s -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out %s"
                  % (input_file, mask_output_file))
        #creating database
        os.system("makeblastdb -in %s -input_type fasta -dbtype nucl -parse_seqids -mask_data %s -out %s -title '%s'"
                  % (input_file, mask_output_file, db_name, db_name))
        #cheking dqatabase
        os.system("blastdbcmd -db %s -info" % db_name)


class BLAST():
    def make_blast_db(self):
        #TODO: write this function
        pass

if __name__ == "__main__":
    """
    #os.chdir("/home/mahajrod/genetics/tetraodon/reference_cda/tetraodon_base")
    BLAST_plus.make_blast_plus_db("/home/mahajrod/genetics/tetraodon/reference_cda/tetraodon_base/assembly_v7.fa.masked",
                  "/home/mahajrod/genetics/tetraodon/reference_cda/tetraodon_base/assembly_v7.asnb",
                  "assembly_v7")
    """
    #os.chdir("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference")
    os.chdir("/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/")
    blast = BLAST_plus()
    """
    blast.make_blast_plus_db("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24s.fasta",
                            "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24s.asnb",
                            "Alsu24s")
    """
    """
    fasta_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.fasta"
    asnb_file = "/home/mahajrod/genetics/desaminases/data/LAN210_v0.10m/LAN210_v0.10m.asnb"
    name = "LAN210_v0.10m"
    """
    os.chdir("/home/mahajrod/genetics/desaminases/data/S288C_R64/")
    fasta_file = "/home/mahajrod/genetics/desaminases/data/S288C_R64/S288C_reference_sequence_R64-1-1_20110203_modified.fsa"
    asnb_file = "/home/mahajrod/genetics/desaminases/data/S288C_R64/S288C_reference_sequence_R64-1-1_20110203_modified.asnb"

    os.chdir("/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/")
    fasta_file = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/Alsu24mc.fasta"
    asnb_file = "/media/mahajrod/d9e6e5ee-1bf7-4dba-934e-3f898d9611c8/Data/Alsu/reference/Alsu24mc/Alsu24mc.asnb"


    name = "Alsu24mc"
    blast.make_blast_plus_db(fasta_file,
                            asnb_file,
                            name)