__author__ = 'mahajrod'

from setuptools import setup, find_packages
from os.path import join, dirname


setup(name='MAVR',
      version='1.0',
      packages=find_packages(),
      author='Sergei F. Kliver',
      author_email='mahajrod@gmail.com',
      long_description=open(join(dirname(__file__), 'README')).read(),)