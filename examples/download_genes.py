import sys
import os
current_dir = os.getcwd()
sys.path.append("../")
from Download.NCBI import *


download_genes_from_geneslist_file("mahajrod@gmail.com", "gene_list.tsv", "genes.gb", "gb")