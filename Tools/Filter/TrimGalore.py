__author__ = 'mahajrod'

from Tools.Abstract import Tool


class TrimGalore(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "trim_galore", path=path, max_threads=max_threads)

    def filter(self, min_length,
               forward_reads,
               forward_trim,
               reverse_reads=None,
               reverse_trim=None,
               quality_score="phred33",
               adapter="AGATCGGAAGAGC",
               adapter2=None,
               quality_treshold=20,
               output_folder="trimmed",
               stringency=None,
               compess=False):
        #TODO: check and add rest of options

        options = ""
        options += " --paired" if reverse_reads else ""
        options += " -a %s" % adapter
        options += " -a2 %s" % adapter2 if adapter2 is not None else ""
        options += " -s %i" % stringency if stringency is not None else ""
        options += " --%s" % quality_score
        options += " --length %i" % min_length
        options += " --%s" % quality_score
        options += " --%s" % quality_score
        options += "" if compess else "--dont_gzip"
        options += " -q %i" % quality_treshold
        options += " -o %s" % output_folder

        options += " --clip_R1 %i" % forward_trim if forward_trim else ""
        options += " --clip_R2 %i" % reverse_trim if reverse_trim else ""

        options += " %s" % forward_reads
        options += " %s" % reverse_reads if reverse_reads else ""

        self.execute(options)