__author__ = 'mahajrod'

from Tools.Annotation.VEP import VEP
from Tools.Annotation.SNPeff import SNPeff
from Tools.Annotation.Barrnap import Barrnap
from Tools.Annotation.Emapper import Emapper
from Tools.Annotation.AUGUSTUS import AUGUSTUS
from Tools.Annotation.Exonerate import Exonerate
from Tools.Annotation.TransDecoder import TransDecoder

SNPeff_path = ""
SNPeff = SNPeff(jar_path=SNPeff_path)

VEP = VEP()
Barrnap = Barrnap()
Emapper = Emapper()
AUGUSTUS = AUGUSTUS()
Exonerate = Exonerate()
TransDecoder = TransDecoder()
