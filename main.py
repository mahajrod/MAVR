#!/usr/bin/env python2

from General.General import *
from MutAnalysis.Mutation import *
#from SeqAnalysis.SeqAnalysis	import *

scr_arg = ScriptArg()

source_files = SourceFiles(scr_arg.arg.source, verbose=scr_arg.arg.verbose)

for filename in source_files.files_dict["vcf"]:
    mutations_vcf = MutationsVcf(filename, from_file=True)
    filtered_mutations, filtered_out_mutations = mutations_vcf.filter_by_reference_and_alt([("G", ["A"]), ("C", ["T"])])
    filtered_mutations.write(filename[:-4] + "_filtered.vcf")
    filtered_out_mutations.write(filename[:-4] + "_filtered_out.vcf")
    #mutations_vcf.write_pos_data()
    #print(mutations_vcf.metadata.metadata)