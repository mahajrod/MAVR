__author__ = 'mahajrod'

from Tools.GATK.ApplyRecalibration import ApplyRecalibration
from Tools.GATK.CombineVariants import CombineVariants
from Tools.GATK.FastaAlternateReferenceMaker import FastaAlternateReferenceMaker
from Tools.GATK.SelectVariants import SelectVariants
from Tools.GATK.UnifiedGenotyper import UnifiedGenotyper
from Tools.GATK.VariantFiltration import VariantFiltration
from Tools.GATK.VariantRecalibrator import VariantRecalibrator

ApplyRecalibration = ApplyRecalibration()
CombineVariants = CombineVariants()
FastaAlternateReferenceMaker = FastaAlternateReferenceMaker()
SelectVariants = SelectVariants()
UnifiedGenotyper = UnifiedGenotyper()
VariantFiltration = VariantFiltration()
VariantRecalibrator = VariantRecalibrator()
