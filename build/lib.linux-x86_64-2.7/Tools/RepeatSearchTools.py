__author__ = 'mahajrod'
import os
from General import check_path
#from Converters.Converters import gff22gff3, gff32gtf


def RepeatModeler_search(query_file, db_name, output_file="run.out",
                         num_of_threads=5, RepeatModeler_dir=""):
    print("\nRepeatModeler search...\n")
    repmod_dir = check_path(RepeatModeler_dir)
    os.system(repmod_dir + "BuildDatabase -engine ncbi  -name %s %s" % (db_name, query_file))
    os.system(repmod_dir + "RepeatModeler -engine ncbi -pa %i -database %s > %s"
              % (num_of_threads, db_name, output_file))


def TRF_search(query_file, match=2, mismatch=7, delta=7, PM=80,
               PI=10, minscore=50, max_period=500, flanked=False, TRF_dir=""):

    print("\nTRF search...\n")
    #use: trf File Match Mismatch Delta PM PI Minscore MaxPeriod [options]
    #Where: (all weights, penalties, and scores are positive)
    # File = sequences input file
    # Match = matching weight
    # Mismatch = mismatching penalty
    # Delta = indel penalty
    # PM = match probability (whole number)
    # PI = indel probability (whole number)
    # Minscore = minimum alignment score to report
    # MaxPeriod = maximum period size to report
    # [options] = one or more of the following :
    # -m masked sequence file
    # -f flanking sequence
    # -d data file
    # -h suppress HTML output
    #Recomended options: trf yoursequence.txt 2 7 7 80 10 50 500 -f -d -m
    flanking = ""
    if flanked:
        flanking = "-f"

    trf_path = check_path(TRF_dir)
    os.system(trf_path + "trf %s %i %i %i %i %i %i %i %s -d -m"
              % (query_file, match, mismatch, delta, PM, PI, minscore, max_period, flanking))


def extract_repbase(species, output_file="RepBase.fasta", RepeatMaskerUtils_dir=""):
    print("\nExtracting RepBase for %s\n" % species)
    repmaskutils_dir = check_path(RepeatMaskerUtils_dir)
    os.system(repmaskutils_dir + "queryRepeatDatabase.pl -species %s > %s" % (species, output_file))


def RepeatMasker_search(query_file, species, custom_lib_path=None, RepeatMasker_dir="",
                        num_of_threads=5, search_type="-s"):

    #species: see list of possible species in repeatmasker.help coming with RepeatMasker
    #search type: "-s" (sensetive), "" (default), "-q" (fast), "-qq" (very fast)

    repmask_dir = check_path(RepeatMasker_dir)
    custom_lib = ""
    if custom_lib_path:
        cuatom_lib = "-lib %s" % custom_lib_path

    #additional options:
    #-xm    creates an additional output file in cross_match format (for parsing)
    #-ace   creates an additional output file in ACeDB format
    #-gff   creates an additional Gene Feature Finding format
    #-excln The percentages displayed in the .tbl file are calculated using a
    #       total sequence length excluding runs of 25 Ns or more.
    print("\nRepeatMasker search...\n")
    os.system(repmask_dir + "RepeatMasker -excln -xm -ace -gff %s -pa %i -species %s %s %s"
              % (custom_lib, num_of_threads, species, search_type, query_file))


def rmout2gff3(rmoutfile, outfile, RepeatMaskerUtils_dir=""):
    repmaskutils_dir = check_path(RepeatMaskerUtils_dir)
    os.system(repmaskutils_dir + "rmOutToGFF3.pl %s > %s" % (rmoutfile, outfile))


def windowmasker_search(windowmasker_dir):
    winmask_dir = check_path(windowmasker_dir)
    #TODO: write this function
    pass

if __name__ == "__main__":
    reference_name = "LAN210_v0.3m"
    #reference_name = "LAN210_v0.4m"
    #reference_name = "LAN210_v0.5m"
    reference_file = reference_name + ".fasta"


    workdir = "/home/mahajrod/genetics/desaminases/data/%s/masking" % reference_name
    os.system("mkdir -p %s" % workdir)
    os.chdir(workdir)
    """
    os.system("ln -fs ../%s %s" % (reference_file, reference_file))
    os.system("mkdir -p repeatmodeler")
    os.chdir("repeatmodeler")
    RepeatModeler_search("../%s" % reference_file, reference_name,
                         RepeatModeler_dir="/home/mahajrod/Repositories/genetic/NGS_tools/RepeatModeler")
    os.chdir(workdir)
    os.system("mkdir -p custom_lib")
    #os.system("/bin/cp -rf ")
    extract_repbase("fungi", output_file="custom_lib/RepBase_fungi.fasta",
                    RepeatMaskerUtils_dir="/home/mahajrod/Repositories/genetic/NGS_tools/RepeatMasker/util")
    """
    #noncontinuos pipeline because of repeatmodeler - it creates diretory with data in name
    #put output of repeatmodeler (consensi.fa.classified) to custom_lib dir and combine it with repbase into custom_lib.fasta
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
