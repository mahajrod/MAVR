__author__ = 'mahajrod'


def check_path(path_to_check):
    #returns path with / at end or blank path
    if path_to_check != "":
        if path_to_check[-1] != "/":
            return path_to_check + "/"
    return path_to_check
