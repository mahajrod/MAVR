
from Pipelines.SNPCall import SNPCallPipeline
from Pipelines.Primers import STRPrimerPipeline
from Pipelines.Abstract import Pipeline
from Pipelines.Filtering import FilteringPipeline
from Pipelines.DiffExpression import DiffExpressionPipeline

Pipeline = Pipeline()
SNPCallPipeline = SNPCallPipeline()
FilteringPipeline = FilteringPipeline()
STRPrimerPipeline = STRPrimerPipeline()
DiffExpressionPipeline = DiffExpressionPipeline()
