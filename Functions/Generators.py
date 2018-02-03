__author__ = 'mahajrod'

from collections import Iterable, OrderedDict


def recursive_generator(iterable):
    if isinstance(iterable, Iterable) and (not isinstance(iterable, str)):
        if isinstance(iterable, dict):
            for element in iterable:
                for ele in recursive_generator(iterable[element]):
                    yield ele
        else:
            for element in iterable:
                for ele in recursive_generator(element):
                    yield ele
    else:
        yield iterable


def recursive_generator_by_type(iterable, type_list):
    if isinstance(iterable, Iterable) and (not isinstance(iterable, str)):
        if isinstance(iterable, dict):
            for element in iterable:
                for ele in recursive_generator(iterable[element]):
                    yield ele
        else:
            for element in iterable:
                for ele in recursive_generator(element):
                    yield ele
    else:
        if type(iterable) in type_list:
            yield iterable
