__author__ = 'mahajrod'


from Tools.Abstract import Tool


class Tophat(Tool):
    def __init__(self, path="", max_threads=4):
        super(Tool).__init__("tophat", path=path, max_threads=max_threads)