from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy

setup(
    name = 'rational',
    ext_modules=[ 
        Extension("rational",  ["rational.pyx"],
        include_dirs = [numpy.get_include()])],
    cmdclass = {'build_ext': build_ext}
)
