#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

import os

from Tools.Abstract import Tool


class Supernova(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "supernova", path=path, max_threads=max_threads)