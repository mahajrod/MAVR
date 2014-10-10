__author__ = 'mahajrod'


from Tools.Abstract import Tool


class Tophat(Tool):
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "tophat", path=path, max_threads=max_threads)