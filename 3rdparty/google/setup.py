#!/usr/bin/env python

import os
import sys
from distutils.core import setup, Extension

os_macro = 'OS_MACOSX' if sys.platform == 'darwin' else 'OS_LINUX'

extension = Extension(
    name='s2',
    sources=['base/int128.cc',
             'base/logging.cc',
             'base/stringprintf.cc',
             'base/strtoint.cc',
             'strings/ascii_ctype.cc',
             'strings/stringprintf.cc',
             'strings/strutil.cc',
             'strings/split.cc',
             'util/coding/coder.cc',
             'util/coding/varint.cc',
             'util/math/mathutil.cc',
             'util/math/mathlimits.cc',
             'util/math/exactfloat/exactfloat.cc',
             's2/s1angle.cc',
             's2/s2.cc',
             's2/s2cellid.cc',
             's2/s2latlng.cc',
             's2/s1interval.cc',
             's2/s2cap.cc',
             's2/s2cell.cc',
             's2/s2cellunion.cc',
             's2/s2edgeindex.cc',
             's2/s2edgeutil.cc',
             's2/s2latlngrect.cc',
             's2/s2loop.cc',
             's2/s2pointregion.cc',
             's2/s2polygon.cc',
             's2/s2polygonbuilder.cc',
             's2/s2polyline.cc',
             's2/s2r2rect.cc',
             's2/s2region.cc',
             's2/s2regioncoverer.cc',
             's2/s2regionintersection.cc',
             's2/s2regionunion.cc',
             's2/s2module.cc',
             's2/s2object.cc',
             's2/s2helper.cc'],
    include_dirs=['.'],
    libraries=['crypto'],
    extra_compile_args=['-std=c++11', '-fno-wrapv'],
    define_macros=[('NDEBUG', '1'),
                   (os_macro, None),
                   ('HASH_NAMESPACE', '__gnu_cxx'),
                   ('ARCH_K8', '1'),
                   ('S2_USE_EXACTFLOAT', None)])

setup(
    name="s2",
    version='0.1',
    description='Google\'s s2 geometry library',
    ext_modules=[extension])
