__author__ = 'mahajrod'

from Tools.Alignment.BWA import BWA
from Tools.Alignment.BLAT import BLAT
from Tools.Alignment.TMAP import TMAP
from Tools.Alignment.STAR import STAR
from Tools.Alignment.Tophat import Tophat
from Tools.Alignment.Bowtie2 import Bowtie2
from Tools.Alignment.Novoalign import Novoalign

max_threads = 4
bowtie2_path = ""
bwa_path = ""
novoalign_path = ""
tmap_path = ""
prank_path = ""


BWA = BWA(path=bwa_path, max_threads=max_threads)
TMAP = TMAP(path=novoalign_path, max_threads=max_threads)
BLAT = BLAT()
STAR = STAR()
Tophat = Tophat()
Bowtie2 = Bowtie2(path=bowtie2_path, max_threads=max_threads)
Novoalign = Novoalign(path=novoalign_path, max_threads=max_threads)
