import os
from multiprocessing import Pool

from Tools.BLAST import *

#TODO: write and rewrite!!!
def transfer_annotations_splign(work_dir,
                                cdna_file,
                                genome_file,
                                blast_dir="/home/mahajrod/Repositories/genetic/BLAST/blast-2.2.9-amd64-linux"):
    os.chdir(work_dir)
    os.system("splign -mklds %s" % work_dir)
    os.system("samtools faidx %s" % cdna_file)
    os.system("samtools faidx %s" % genome_file)

    print("Creating BLAST database from cdna sequences")
    os.system(blast_dir + "/formatdb -pF -oT -i %s" % cdna_file)
    print("Creating BLAST database from genome sequences")
    os.system(blast_dir + "/formatdb -pF -oT -i %s" % genome_file)

    #make_blast_db(cdna_file, cdna_file[:-3] + ".asnb", cdna_db_name)

    #make_blast_db(genome_file, genome_file[:-7] + ".asnb", genome_db_name)

    os.system("compart -qdb %s -sdb %s > cdna.compartments" % (cdna_file, genome_file))
    print("Running splign")
    os.system("splign -ldsdir %s -comps cdna.compartments > splign.out" % (work_dir))
    """
    print("Creating BLAST database from cdna sequences")
    os.system(blast_dir + "/formatdb -i %s -p F" % genome_file)
    os.system(blast_dir + "/megablast -i %s -d %s -F 'm D;R' -D 3 | grep -v '^#' | sort -k 2,2 -k 1,1 -T temp_dir > cdna.hit" %
              (cdna_file, genome_file))
    os.system("splign -ldsdir %s -hits cdna.hit > splign.out" % work_dir)
    """


def transfer_annotations_prosplign(work_dir,
                                protein_file,
                                genome_blast_db,
                                blast_dir="/home/mahajrod/Repositories/genetic/BLAST/blast-2.2.9-amd64-linux"):

    os.system("tblastn -query %s -db %s -outfmt 6 | sort -k 2,2 -k 1,1 > tblastn.hit")
    os.system("compart -f blast.hit -add 10000 > comp")
    #prosplign -two_stages -pfa p.fa -nfa n.fa -f comp -inf pro.inf -out pro.out

def arg_splitter(args):
    return transfer_annotations_splign(*args)

if __name__ == "__main__":
    reference_species = [
                        ("astyanax_mexicanus", "Astyanax_mexicanus.AstMex102.75.cdna.all.fa"),
                        ("danio_rerio", "Danio_rerio.Zv9.75.cdna.all.fa"),
                        ("gadus_morhua", "Gadus_morhua.gadMor1.75.cdna.all.fa"),
                        ("gasterosteus_aculeatus", "Gasterosteus_aculeatus.BROADS1.75.cdna.all.fa"),
                        ("latimeria_chalumnae", "Latimeria_chalumnae.LatCha1.75.cdna.all.fa"),
                        ("oreochromis_niloticus", "Oreochromis_niloticus.Orenil1.0.75.cdna.all.fa"),
                        ("oryzias_latipes", "Oryzias_latipes.MEDAKA1.75.cdna.all.fa"),
                        ("takifugu_rubripes", "Takifugu_rubripes.FUGU4.75.cdna.all.fa"),
                        ("xiphophorus_maculatus", "Xiphophorus_maculatus.Xipmac4.4.2.75.cdna.all.fa")
                        ]

    work_dir = "/home/mahajrod/genetics/tetraodon/reference_cda"
    arg_list = []

    for species, cdna_file in reference_species:
        arg_list.append((work_dir + "/" + species + "/cdna",
                         work_dir + "/" + species + "/cdna/" + cdna_file,
                         work_dir + "/" + species + "/cdna/assembly_v7.fa.masked"))
        """
        transfer_annotations_splign(work_dir + "/" + species + "/cdna",
                                    work_dir + "/" + species + "/cdna/" + cdna_file,
                                    work_dir + "/" + species + "/cdna/assembly_v7.fa.masked")
        """
    pool = Pool(processes=4)
    pool.map(arg_splitter, arg_list)
    """
    transfer_annotations_splign("/home/mahajrod/genetics/tetraodon/reference_cda/astyanax_mexicanus/cdna",
                                "/home/mahajrod/genetics/tetraodon/reference_cda/astyanax_mexicanus/cdna/Astyanax_mexicanus.AstMex102.75.cdna.all.fa",
                                "AstMex102.75.cdna",
                                "/home/mahajrod/genetics/tetraodon/reference_cda/astyanax_mexicanus/cdna/assembly_v7.fa.masked",
                                "assembly_v7")
    """
