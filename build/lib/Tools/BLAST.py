__author__ = 'mahajrod'
import os


def make_blast_plus_db(input_file, mask_output_file, db_name):
    #makes BLAST database from fasta file
    os.system("dustmasker -in %s -infmt fasta -parse_seqids -outfmt maskinfo_asn1_bin -out %s"
              % (input_file, mask_output_file))
    #creating database
    os.system("makeblastdb -in %s -input_type fasta -dbtype nucl -parse_seqids -mask_data %s -out %s -title '%s'"
              % (input_file, mask_output_file, db_name, db_name))
    #cheking dqatabase
    os.system("blastdbcmd -db %s -info" % db_name)


def make_blast_db():
    #TODO: write this function
    pass

if __name__ == "__main__":
    os.chdir("/home/mahajrod/genetics/tetraodon/reference_cda/tetraodon_base")
    make_blast_plus_db("/home/mahajrod/genetics/tetraodon/reference_cda/tetraodon_base/assembly_v7.fa.masked",
                  "/home/mahajrod/genetics/tetraodon/reference_cda/tetraodon_base/assembly_v7.asnb",
                  "assembly_v7")