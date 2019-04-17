#/usr/bin/env python

from Bio import Entrez
from RouToolPa.Parsers.General import parse_sv

"""
def download_genes_by_name(email, gene_name, output_filename, output_type):
    Entrez.email = email
    search = Entrez.esearch(db="gene", term=gene_name+"[sym]", retmax=10000, rettype="tab")
    handled_search = Entrez.read(search)
    #handle_search = Entrez.efetch(db="nucleotide", id="186972394", rettype="gb", retmode="text")
    with open("search.tsv", "w") as search_fd:
        for gene_id in handled_search["IdList"]:
            handle = Entrez.esummary(db="gene", id=gene_id, retmax=10000, rettype=)
	        #seq_handle = Entrez.efetch(db="gene", id=gene_id, rettype="tab", retmode="txt")




        #record = Entrez.read(Entrez.elink(dbfrom="gene", id=gene_id))
        record = Entrez.read(Entrez.elink(db="nucleotide", dbfrom="gene", id=gene_id))
        print(gene_id, record)#[0]["LinkSetDb"][0]["Link"][0])
        #record = Entrez.read(ret_search)
    #print(record.keys())
    #print(record)
"""


def download_genes_from_geneslist_file(email, input_file, output_file, output_type):
    Entrez.email = email
    header, record_list = parse_sv(input_file)
    print("Downloading %i genes" % len(record_list))
    with open(output_file, "w") as fd:
        pass
    i = 1
    for record in record_list:

        source_id = record[11]
        start = int(record[12])
        end = int(record[13])
        strand = record[14]
        seq_handle = Entrez.efetch(db="nucleotide", id=source_id, seq_start=start, seq_stop=end, rettype=output_type, retmode="txt")
        with open(output_file, "a") as fd:
            fd.write(seq_handle.read())
        print("Downloaded %i from %i genes" % (i, len(record_list)))
        i += 1