__author__ = 'mahajrod'

from Tools.Filter.FaCut import FaCut
from Tools.Filter.FastQC import FastQC
from Tools.Filter.Adapters import adapters_PE
from Tools.Filter.Cutadapt import Cutadapt
from Tools.Filter.Clumpify import Clumpify
from Tools.Filter.TrimGalore import TrimGalore
from Tools.Filter.Trimmomatic import Trimmomatic
from Tools.Filter.Cookiecutter import Cookiecutter
from Tools.Filter.Cookiecutter import CookiecutterOld

FaCut = FaCut()
FastQC = FastQC()
Cutadapt = Cutadapt()
Clumpify = Clumpify()
TrimGalore = TrimGalore()
Trimmomatic = Trimmomatic()
Cookiecutter = Cookiecutter()
CookiecutterOld = CookiecutterOld()

