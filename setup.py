from setuptools import setup, find_packages, Command
from setuptools.command.develop import develop
from setuptools.command.install import install
from setuptools.command.egg_info import egg_info
import atexit, os, sys

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
import shutil
import numpy as np
from glob import glob

ext_modules = [
    Extension("Delta", ["SGVFinder2/helpers/Delta.pyx"]),
    Extension("sam2pmp_helper", ["SGVFinder2/helpers/sam2pmp_helper.pyx"]),
]

setup(
    author="Tal Korem",
    author_email="tk2829@cumc.columbia.edu",
    version="v0.0.1",
    name="SGVFinder2",
    package_dir={"SGVFinder2": "SGVFinder2"},
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "numpy",
        "pandas",
        "ujson",
        "cython",
        "pysam",
        "scipy",
        "bokeh"
    ],
    python_requires=">=3",
    license="Apache-2.0",
    license_files=["LICENSE"],
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules
)
