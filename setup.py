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
    version="1.1.0",
    url="https://www.github.com/pirovc/grimer",
    license='MIT',
    author="Vitor C. Piro",
    author_email="pirovc@posteo.net",
    description="GRIMER: contamination detection and microbiome exploration",
    long_description=read("README.md"),
    packages=['grimer'],
    package_data={
        'grimer': ['css/*', 'img/*', 'js/*', 'scripts/*']
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
