
from Tools.Abstract import Tool


class Cufflinks(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "cufflinks", path=path, max_threads=max_threads)