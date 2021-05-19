# Copyright (c) 2017-2021 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

import os
import re
import sys
import setuptools
from setuptools import setup, Extension
from numpy import get_include

sys.path.append(os.path.join(os.path.dirname(__file__), 'repack'))
from VERSION import __version__


srcdir = './repack/utils/'  # C-code source folder

files = os.listdir(srcdir)
files = list(filter(lambda x: re.search('.+[.]c$', x), files))
files = list(filter(lambda x: not re.search('[.#].+[.]c$', x), files))

inc = [get_include()]
eca = ['-ffast-math']
ela = []

extensions = [
    Extension(
        'repack.utils.' + cfile.rstrip('.c'),
        sources=[f'{srcdir}{cfile}'],
        include_dirs=inc,
        extra_compile_args=eca,
        extra_link_args=ela,
        )
    for cfile in files
    ]

with open('README.md', 'r') as f:
    readme = f.read()


setup(
    name = 'lbl-repack',
    version = __version__,
    author = 'Patricio Cubillos',
    author_email = 'patricio.cubillos@oeaw.ac.at',
    url = 'https://github.com/pcubillos/repack',
    packages = setuptools.find_packages(),
    install_requires = [
        'numpy>=1.13.3',
        'scipy>=0.17.1',
        ],
    include_package_data = True,
    license = 'MIT',
    description = 'A line-transition data compression package.',
    long_description = readme,
    long_description_content_type = 'text/markdown',
    include_dirs = inc,
    entry_points = {'console_scripts': ['repack = repack.__main__:main']},
    ext_modules = extensions,
    )
