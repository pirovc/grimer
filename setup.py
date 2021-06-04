import io
import os
import re

from setuptools import setup


def read(filename):
    filename = os.path.join(os.path.dirname(__file__), filename)
    text_type = type(u"")
    with io.open(filename, mode="r", encoding='utf-8') as fd:
        return re.sub(text_type(r':[a-z]+:`~?(.*?)`'), text_type(r'``\1``'), fd.read())

setup(
    name="grimer",
    version="0.0.0",
    url="https://www.github.com/pirovc/grimer",
    license='MIT',

    author="Vitor C. Piro",
    author_email="pirovc@posteo.net",

    description="GRIMER portable microbiome visualization and contamination detection",
    long_description=read("README.md"),

    packages=['grimer'],
    #install_requires=['binpacking==1.4.3'],

    include_package_data=True,
    package_data={
        'js': ['*'],
        'css': ['*'],
        'img': ['*']
    },

    entry_points={'console_scripts': ['grimer=grimer.grimer:main']},

    classifiers=[
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
