__author__ = 'mahajrod'

from Tools.GATK.IndelRealigner import IndelRealigner
from Tools.GATK.SelectVariants import SelectVariants
from Tools.GATK.CombineVariants import CombineVariants
from Tools.GATK.BaseRecalibrator import BaseRecalibrator
from Tools.GATK.UnifiedGenotyper import UnifiedGenotyper
from Tools.GATK.VariantFiltration import VariantFiltration
from Tools.GATK.ApplyRecalibration import ApplyRecalibration
from Tools.GATK.VariantRecalibrator import VariantRecalibrator
from Tools.GATK.RealignerTargetCreator import RealignerTargetCreator
from Tools.GATK.FastaAlternateReferenceMaker import FastaAlternateReferenceMaker


IndelRealigner = IndelRealigner()
SelectVariants = SelectVariants()
CombineVariants = CombineVariants()
BaseRecalibrator = BaseRecalibrator()
UnifiedGenotyper = UnifiedGenotyper()
VariantFiltration = VariantFiltration()
ApplyRecalibration = ApplyRecalibration()
VariantRecalibrator = VariantRecalibrator()
RealignerTargetCreator = RealignerTargetCreator()
FastaAlternateReferenceMaker = FastaAlternateReferenceMaker()
