from numpy import get_include
import os, re, sys
from setuptools import setup, Extension

topdir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(topdir + "/repack")
import VERSION as ver

srcdir = './repack/utils/'  # C-code source folder

# Get all file from source dir:
files = os.listdir(srcdir)

# This will filter the results for just the c files:
files = filter(lambda x:     re.search('.+[.]c$',     x), files)
files = filter(lambda x: not re.search('[.#].+[.]c$', x), files)

inc = [get_include()]
eca = ['-ffast-math']
ela = []

extensions = []
for i in range(len(files)):
  e = Extension(files[i].rstrip('.c'),
                sources=["{:s}{:s}".format(srcdir, files[i])],
                include_dirs=inc,
                extra_compile_args=eca,
                extra_link_args=ela)
  extensions.append(e)

setup(name         = "repack",
      version      = '{:d}.{:d}.{:d}'.format(ver.repack_VER,
                                             ver.repack_MIN, ver.repack_REV),
      author       = "Patricio Cubillos",
      author_email = "patricio.cubillos@oeaw.ac.at",
      url          = "https://github.com/pcubillos/repack",
      #packages     = ["repack"],
      license      = ["MIT"],
      description  = 'Re-pack line-transition data.',
      include_dirs = inc,
      ext_modules  = extensions)

