from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
  name = 'nanosleep',
  ext_modules=[ 
    Extension("nanosleep",  ["nanosleep.pyx"])],
    cmdclass = {'build_ext': build_ext}
)
