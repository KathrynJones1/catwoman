import setuptools
from setuptools import setup, Extension, find_packages
import numpy as np
import re
import os
import sys
import tempfile
from distutils.ccompiler import new_compiler

# Detect OpenMP support (using distutils.ccompiler is fine here)
def hasfunction(cc, funcname, include=None, extra_postargs=None):
    tmpdir = tempfile.mkdtemp(prefix='hasfunction-')
    devnull = oldstderr = None
    try:
        try:
            fname = os.path.join(tmpdir, 'funcname.c')
            with open(fname, 'w') as f:
                if include is not None:
                    f.write(f'#include {include}\n')
                f.write('int main(void) {\n')
                f.write(f'    {funcname};\n')
                f.write('}\n')
            devnull = open(os.devnull, 'w')
            oldstderr = os.dup(sys.stderr.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            objects = cc.compile([fname], output_dir=tmpdir, extra_postargs=extra_postargs)
            cc.link_executable(objects, os.path.join(tmpdir, "a.out"))
        except Exception:
            return False
        return True
    finally:
        if oldstderr is not None:
            os.dup2(oldstderr, sys.stderr.fileno())
        if devnull is not None:
            devnull.close()

def detect_openmp():
    compiler = new_compiler()
    print("Checking for OpenMP support... ")
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads()')
    needs_gomp = hasopenmp
    if not hasopenmp:
        compiler.add_library('gomp')
    hasopenmp = hasfunction(compiler, 'omp_get_num_threads()')
    needs_gomp = hasopenmp
    if hasopenmp:
        print("Compiler supports OpenMP")
    else:
        print("Did not detect OpenMP support.")
    return hasopenmp, needs_gomp

has_openmp, needs_gomp = detect_openmp()
parallel_args = ['-fopenmp', '-std=c99'] if has_openmp else ['-std=c99']
parallel_libraries = ['gomp'] if needs_gomp else []

ext_modules = [
    Extension('catwoman._nonlinear_ld', ['c_src/_nonlinear_ld.c'],
              include_dirs=[np.get_include()],
              extra_compile_args=parallel_args, libraries=parallel_libraries),
    Extension('catwoman._quadratic_ld', ['c_src/_quadratic_ld.c'],
              include_dirs=[np.get_include()],
              extra_compile_args=parallel_args, libraries=parallel_libraries),
    Extension('catwoman._logarithmic_ld', ['c_src/_logarithmic_ld.c'],
              include_dirs=[np.get_include()],
              extra_compile_args=parallel_args, libraries=parallel_libraries),
    Extension('catwoman._exponential_ld', ['c_src/_exponential_ld.c'],
              include_dirs=[np.get_include()],
              extra_compile_args=parallel_args, libraries=parallel_libraries),
    Extension('catwoman._custom_ld', ['c_src/_custom_ld.c'],
              include_dirs=[np.get_include()],
              extra_compile_args=parallel_args, libraries=parallel_libraries),
    Extension('catwoman._power2_ld', ['c_src/_power2_ld.c'],
              include_dirs=[np.get_include()],
              extra_compile_args=parallel_args, libraries=parallel_libraries),
    Extension('catwoman._rsky', ['c_src/_rsky.c'],
              include_dirs=[np.get_include()],
              extra_compile_args=parallel_args, libraries=parallel_libraries),
    Extension('catwoman._eclipse', ['c_src/_eclipse.c'],
              include_dirs=[np.get_include()],
              extra_compile_args=parallel_args, libraries=parallel_libraries),
]

VERSIONFILE='catwoman/__init__.py'
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo: 
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(
    name='catwoman',
    version=verstr,
    author='Kathryn Jones',
    author_email='kathryndjones@hotmail.co.uk',
    url='https://github.com/KathrynJones1/catwoman',
    packages=find_packages(),
    license='GNU GPLv3',
    description='Transit modelling package for asymmetric light curves',
    long_description=open('README.rst', encoding='utf-8').read(),
    long_description_content_type='text/x-rst',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering',
        'Programming Language :: Python',
    ],
    install_requires=['numpy>=1.16.2'],
    setup_requires=['wheel', 'numpy>=1.16.2'],
    extras_require={
        'matplotlib': ['matplotlib'],
    },
    ext_modules=ext_modules,
)
