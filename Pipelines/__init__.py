from Pipelines.Sanger import SangerPipeline
from Pipelines.SNPCall import SNPCallPipeline
from Pipelines.Primers import STRPrimerPipeline
from Pipelines.Primers import MitochondrialAmplificationPrimerPipeline
from Pipelines.Abstract import Pipeline
from Pipelines.Alignment import AlignmentPipeline
from Pipelines.Filtering import FilteringPipeline
from Pipelines.DiffExpression import DiffExpressionPipeline

Pipeline = Pipeline()
SangerPipeline = SangerPipeline()
SNPCallPipeline = SNPCallPipeline()
AlignmentPipeline = AlignmentPipeline()
FilteringPipeline = FilteringPipeline()
STRPrimerPipeline = STRPrimerPipeline()
DiffExpressionPipeline = DiffExpressionPipeline()
MitochondrialAmplificationPrimerPipeline = MitochondrialAmplificationPrimerPipeline()
