from distutils.core import setup
from distutils.extension import Extension

import os
import numpy

try:
    from Cython.Build import cythonize
    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

if not os.path.exists('grizli/utils_c/interp.pyx'):
    USE_CYTHON = False
    
if USE_CYTHON:
    cext = '.pyx'
else:
    cext = '.c'

print('C extension: {0}'.format(cext))

extensions = []
#     Extension("grizli/utils_c/interp", ["grizli/utils_c/interp"+cext],
#         include_dirs = [numpy.get_include()],),
#         
#     # Extension("grizli/utils_c/nmf", ["grizli/utils_c/nmf"+cext],
#     #     include_dirs = [numpy.get_include()],),
#     
#     Extension("grizli/utils_c/disperse", ["grizli/utils_c/disperse"+cext],
#         include_dirs = [numpy.get_include()],),
# 
# ]

if USE_CYTHON:
    extensions = cythonize(extensions)

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "sgas",
    version = "0.0.0",
    author = "Gabriel Brammer",
    author_email = "gbrammer@gmail.com",
    description = "Analysis of the SGAS lens",
    license = "MIT",
    url = "https://github.com/gbrammer/sgas",
    download_url = "https://github.com/gbrammer/sgas/tarball/0.0.0",
    packages=['sgas'], #, 'grizli/utils_c'],
    # requires=['numpy', 'scipy', 'astropy', 'drizzlepac', 'stwcs'],
    # long_description=read('README.rst'),
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Astronomy',
    ],
    ext_modules = extensions,
    package_data={'sgas': ['data/*']},
    # scripts=['grizli/scripts/flt_info.sh'],
)
