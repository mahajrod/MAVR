__author__ = 'mahajrod'

from Tools.Filter.FaCut import FaCut
from Tools.Filter.FastQC import FastQC
from Tools.Filter.Adapters import adapters_PE
from Tools.Filter.TrimGalore import TrimGalore
from Tools.Filter.Trimmomatic import Trimmomatic
from Tools.Filter.Cookiecutter import Cookiecutter


FaCut = FaCut()
FastQC = FastQC()
TrimGalore = TrimGalore()
Trimmomatic = Trimmomatic()
Cookiecutter = Cookiecutter()

