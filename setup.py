#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for pyDataCube
"""

import sys
from setuptools import setup, find_packages
#import glob
import versioneer

#scripts = glob.glob("scripts/*")

setup(name='datacube',
      packages=find_packages(),
      version=versioneer.get_version(),
      mdclass=versioneer.get_cmdclass(),
      description="Cope with the sub-mm observations by using fits files",
      author="Gao-Yuan Zhang",
      author_email="zgy0106@gmail.com",
      setup_requires=["numpy","scipy","matplotlib"],
      install_requires=["numpy","astropy>=2.0","scipy"],
      include_package_data=True,
#      scripts=scripts.
      classifiers=[
          'Intended Audience :: Science/Research',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Visualization',
      ],
)
