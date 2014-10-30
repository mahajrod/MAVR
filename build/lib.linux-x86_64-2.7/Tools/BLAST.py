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


class BLAST(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "blastall", path=path, max_threads=max_threads)

    def make_db(self, input_files, title, input_type="nucleotide", logfile=None, parse_seqid=None, seq_entry_input=None,
                blast_base_names=None):
        """
      -t  Title for database file [String]  Optional
      -i  Input file(s) for formatting [File In]  Optional
      -l  Logfile name: [File Out]  Optional
        default = formatdb.log
      -p  Type of file
             T - protein
             F - nucleotide [T/F]  Optional
        default = T
      -o  Parse options
             T - True: Parse SeqId and create indexes.
             F - False: Do not parse SeqId. Do not create indexes.
     [T/F]  Optional
        default = F
      -a  Input file is database in ASN.1 format (otherwise FASTA is expected)
             T - True,
             F - False.
     [T/F]  Optional
        default = F
      -b  ASN.1 database in binary mode
             T - binary,
             F - text mode.
     [T/F]  Optional
        default = F
      -e  Input is a Seq-entry [T/F]  Optional
        default = F
      -n  Base name for BLAST files [String]  Optional
      -v  Database volume size in millions of letters [Integer]  Optional
        default = 4000
      -s  Create indexes limited only to accessions - sparse [T/F]  Optional
        default = F
      -V  Verbose: check for non-unique string ids in the database [T/F]  Optional
        default = F
      -L  Create an alias file with this name
            use the gifile arg (below) if set to calculate db size
            use the BLAST db specified with -i (above) [File Out]  Optional
      -F  Gifile (file containing list of gi's) [File In]  Optional
      -B  Binary Gifile produced from the Gifile specified above [File Out]  Optional
      -T  Taxid file to set the taxonomy ids in ASN.1 deflines [File In]  Optional
        """
        options = ""
        options += " -t %s" % title
        options += " -l %s" % logfile if logfile else ""
        options += " -p %s" % ("T" if input_type == "protein" else "F")
        options += " -o T" if parse_seqid else ""
        options += " -e T" if seq_entry_input else ""
        options += " -n %s" % blast_base_names if blast_base_names else ""
        # checking if inputfiles has method __iter__
        options += " -i %s" % (" ".join(input_files) if hasattr(input_files, '__iter__') and not isinstance(input_files, str) else input_files)

        #TODO: add parsing of rest of arguments
        self.execute(options, cmd="formatdb")


if __name__ == "__main__":

    blast_plus = BLAST_plus()
    workdir = "/home/mahajrod/genetics/nxf/annotation/test/Dmel/nxf1/"
    os.chdir(workdir)
    blast_plus.make_blast_plus_db("/home/mahajrod/genetics/nxf/annotation/test/Dmel/nxf1/dmel-all-chromosome-r5.54.fasta",
                                  "dmel_all_chromosome.asnb", "dmel_all_chromosome")
    """
    workdir = "/home/mahajrod/genetics/reference/"
    blast = BLAST()

    os.chdir(workdir)
    files = os.listdir("./")

    links_directory = "/home/mahajrod/genetics/nxf/annotation/test/Dmel/nxf1/"
    drosophila_dirs = []
    for entry in files:
        if "drosophila" in entry:
            drosophila_dirs.append(entry)

    for directory in drosophila_dirs:
        os.chdir(workdir)
        os.chdir(directory)
        reference = os.listdir("./")[0] + "/"
        os.chdir(reference)
        os.chdir("ref")
        for file_entry in os.listdir("./"):
            if file_entry[-6:] == ".fasta":
                reference_file = file_entry
        blast.make_db(reference_file, directory, parse_seqid=True)
        #print("ln -s %s %s" % (workdir + directory + "/ref/" + reference_file, links_directory))
        #os.system("ln -fs %s %s" % (workdir + directory + "/" + reference + "ref/" + reference_file, links_directory))
    """