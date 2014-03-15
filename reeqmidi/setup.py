
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    name = 'reeqmidi',
    ext_modules=[ 
        Extension("reeqmidi",  ["reeqmidi.pyx"],
        include_dirs = [numpy.get_include()])],
    cmdclass = {'build_ext': build_ext}
)
