#!/usr/bin/env python
import os
from Tools.Abstract import Tool


class PRANK(Tool):
    def __init__(self, path="", max_threads=1):
        Tool.__init__(self, "prank", path=path, max_threads=max_threads)

    def align(self,
              sequence_file, output, tree_file=None, output_format=None, show_xml=None,
              show_tree=None, show_ancestral_sequences=None, show_evolutionary_events=None,
              showall=None, compute_posterior_support=None, njtree=None, skip_insertions=False,
              codon_alignment=None, translated_alignment=None):
        # TODO: add rest of options
        options = " -d=%s" % sequence_file
        options += " -o=%s" % output

        options += " -t=%s" % tree_file if tree_file else ""
        options += " -f=%s" % output_format if output_format else ""
        options += " -showxml" if not show_xml else ""
        options += " -showtree" if not show_tree else ""
        options += " -showanc" if not show_ancestral_sequences else ""
        options += " -showevents" if not show_evolutionary_events else ""
        options += " -showall" if not showall else ""
        options += " -support" if not compute_posterior_support else ""
        options += " -njtree" if not njtree else ""
        options += " -F" if skip_insertions else ""
        options += " -codon" if codon_alignment else ""
        options += " -translate" if translated_alignment else ""

        self.execute(options)
