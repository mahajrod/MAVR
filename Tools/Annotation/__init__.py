__author__ = 'mahajrod'

from Tools.Annotation.SNPeff import SNPeff
from Tools.Annotation.TransDecoder import TransDecoder
from Tools.Annotation.AUGUSTUS import AUGUSTUS
from Tools.Annotation.Exonerate import Exonerate

SNPeff_path = ""
SNPeff = SNPeff(jar_path=SNPeff_path)

TransDecoder = TransDecoder()
AUGUSTUS = AUGUSTUS()
Exonerate = Exonerate()
