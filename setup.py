#!/usr/bin/env python

import os
import sys
from distutils.core import setup, Extension

os_macro = 'OS_MACOSX' if sys.platform == 'darwin' else 'OS_LINUX'

extension = Extension(
    name='rugis',
    sources=[
        "gpsconv.cc",
        "s2helper.cc",
        "radar.cc",
        "geo.cc",
        "rugis_python.cc"
        ],
    include_dirs=['/usr/lib','/usr/local/lib','3rdparty', '3rdparty/google','3rdparty/boost_1_58_0/include'],
    libraries=['s2', 'crypto'],
    library_dirs=['3rdparty/lib'],
    runtime_library_dirs=['.', './rugis_libs','/usr/lib','/usr/local/lib'],
    extra_compile_args=['-std=c++11', '-fno-wrapv', '-stdlib=libc++'],
    define_macros=[('NDEBUG', '1'), (os_macro, None)]
    )

setup(
    name="rugis",
    version='2.0',
    description='preprocess module for project radar',
    ext_modules=[extension])
