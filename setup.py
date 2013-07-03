#!/usr/bin/env python
from distutils.core import setup, Extension
modulel = Extension('pyclustal',
                    include_dirs = ['clustal-omega/'],
                    libraries = ['clustalo'],
                    library_dirs=['clustal-omega/'],
                    #extra_compile_args = ['-fopenmp'],
                    #extra_link_args = ['-fopenmp'],
                    sources = ['pyclustal.c'])

setup (name = 'pyclustal',
        version = '0.1',
        description = 'This is the python building of clustal-omega',
        author ='Koonwah Chen',
        author_email = 'cadmuxe@gmail',
        long_description = '''
        Clustal:Multiple alignment of nucleic acid and protein sequences.(http://www.clustal.org/)
        ''',
        ext_modules=[modulel])
