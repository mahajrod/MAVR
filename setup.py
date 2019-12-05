__author__ = 'mahajrod'
import os
import sys
from setuptools import setup, find_packages


def get_list_of_files(list_of_dirs_and_files, expression=None, recursive=True):
    file_list = []

    for entry in [list_of_dirs_and_files] if isinstance(list_of_dirs_and_files, str) else list_of_dirs_and_files:

        if os.path.isdir(entry):

            files_in_dir = ["%s%s" % (entry if entry[-1] == "/" else (entry+"/"), filename)
                            for filename in sorted(filter(expression, os.listdir(entry))
                                                   if expression else os.listdir(entry))]
            if recursive:
                for filename in files_in_dir:
                    if os.path.isdir(filename):
                        file_list += get_list_of_files([filename], expression=expression, recursive=True)
                    else:
                        file_list.append(filename)
            else:
                file_list += files_in_dir
        elif os.path.exists(entry):
            file_list.append(os.path.abspath(entry))
        else:
            print("%s does not exist" % entry)
    return file_list


dependencies = ['scipy', 'numpy', 'pandas', 'matplotlib', 'matplotlib-venn',
                'biopython', 'xmltodict', 'bcbio-gff', 'statsmodels']
if sys.version_info[0] == 3:
    dependencies += ["ete3"]
elif sys.version_info[0] == 2:
    dependencies += ["ete2"]
else:
    raise ValueError("ERROR!!! Unsupported python version: %s" % str(sys.version_info[0]))

setup(name='MAVR',
      version='0.93',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      install_requires=dependencies,
      scripts=get_list_of_files("scripts/"),
      long_description=open(os.path.join(os.path.dirname(__file__), 'README.md')).read(),)
