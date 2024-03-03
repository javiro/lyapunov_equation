# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

__VERSION__ = open('VERSION').read().strip()

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name=lyapunov_equation,
    version=__VERSION__,
    description='',
    author='',
    author_email='',
    license='',
    url='',
    long_description='',
    install_requires=required,
    packages=find_packages(exclude=('tests', 'tests.*',)),
    include_package_data=True
)
