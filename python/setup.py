import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from sys import platform

pkg_name = "bie-solvers2dpy"

list_files=['../src/helm_wrappers/pseudo_diff_op2d.f']

FLIBS = os.getenv('FMMBIE_LIBS')
FLIBS = FLIBS.rstrip().split(' ')
FLIBS = list(filter(None,FLIBS))
FLIBS.append('../lib-static/libfmm2dbie.a')
if platform == "darwin":
    FLIBS.append('-L/usr/local/lib')
if platform == "linux" or platform == "linux2":
    FLIBS.append('-lopenblas')

ext_helm = Extension(
    name='fmm2dbie',
    sources=list_files,
    extra_f90_compile_args=["-std=legacy"],
    extra_link_args=FLIBS
)

setup(
    name=pkg_name,
    version="0.1.0",
    # author="Manas Rachh",
    # author_email="mrachh@flatironinstitute.org",
    # description="This pacakge contains basic routines Helmholtz dirichlet fast direct solver",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
        "pytest"
    ],
    ext_modules=[ext_helm],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )    
)