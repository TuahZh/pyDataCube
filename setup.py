#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Setup file for pyDataCube
"""

import sys
from setuptools import setup, find_packages


def setup_package():
    needs_sphinx = {'build_sphinx', 'upload_docs'}.intersection(sys.argv)
    sphinx = ['sphinx'] if needs_sphinx else[]
    setup(#setup_requires=sphinx,
        entry_points=entry_points,
        use_pyscaffold=True)
