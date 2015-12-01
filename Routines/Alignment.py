__author__ = 'mahajrod'

from Bio import SearchIO


class AlignmentRoutines:
    def __init__(self):
        pass

    @staticmethod
    def get_db_ids(search_dict):
        id_set = set()
        for query_id in search_dict:
            for hit in search_dict[query_id]:
                id_set.add(hit.id)
        return id_set




