# -*- coding: utf-8 -*-
__author__ = 'Brice Olivier'

import os
from setuptools import setup


def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()


with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(name='em-analysis',
      version='0.2a0',
      description='A Python tool for analysis eye-movements.',
      install_requires=requirements,
      author=['Brice Olivier', 'Jean-Baptiste Durand'],
      author_email='briceolivier1409@gmail.com',
      url='https://github.com/PyENE/em-analysis/',
      packages=['ema'],
      license=read('LICENSE'),
      long_description=read('README.md'))
