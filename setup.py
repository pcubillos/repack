# Copyright (c) 2017-2024 Patricio Cubillos and contributors.
# repack is open-source software under the MIT license (see LICENSE).

import os
import re
from setuptools import setup, Extension

from numpy import get_include


srcdir = './repack/utils/'  # C-code source folder

files = os.listdir(srcdir)
files = list(filter(lambda x: re.search('.+[.]c$', x), files))
files = list(filter(lambda x: not re.search('[.#].+[.]c$', x), files))

inc = [get_include()]
eca = ['-lm', '-O3', '-ffast-math']
ela = ['-lm']

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

setup(
    ext_modules = extensions,
    include_dirs = inc,
)
