import os
from numpy import f2py
from numpy.distutils.core import Extension, setup

#f_args = ' --fcompiler=intelem'' --f90flags=''-fopenmp'' -lgomp'
f_args = ' --f90flags=-fopenmp -lgomp'

f90_dir = './util/'
f90_files = ('mod_diff.f90', 'mod_metrics.f90', 'interp2d.f90', 'unstructured.f90', 'gridgen.f90', 'flow.f90')

os.chdir(os.path.split(os.path.abspath(__file__))[0] )
os.chdir(f90_dir)

#ext = Extension(name='util_f', sources=f90_files)
#setup(name='util_f', author='Yuchen Liu', ext_modules=[ext])
#exit()

#try:
#    source = b''
#    for filename in f90_files:
#        with open(filename, mode='rb') as fh:
#            source += fh.read()
#    f2py.compile(source, modulename='util_f', extension='.f90', extra_args=f_args)
#except:
#    print('Error in compiling!')
#    exit()


#try:
#    with open("./interp2d.f90") as sourcefile:
#        sourcecode = sourcefile.read()
#        f2py.compile(sourcecode, modulename='interp2d', extension='.f90', extra_args=f_args)
#except:
#    print('Error in compileing interp2d!')
#    exit()
#
#try:
#    with open("./unstructured.f90") as sourcefile:
#        sourcecode = sourcefile.read()
#        f2py.compile(sourcecode, modulename='uns', extension='.f90', extra_args=f_args)
#except:
#    print('Error in compileing uns!')
#    exit()
#
#try:
#    with open("./gridgen.f90") as sourcefile:
#        sourcecode = sourcefile.read()
#        f2py.compile(sourcecode, modulename='gridgen', extension='.f90', extra_args=f_args)
#except:
#    print('Error in compileing gridgen!')
#    exit()
#
##try:
##    with open("./mod_metrics.f90") as sourcefile:
##        sourcecode = sourcefile.read()
##        f2py.compile(sourcecode, modulename='mm', extension='.f90', extra_args=f_args)
##except:
##    print 'Error in compileing modmetrics!'
##    exit()
#
#try:
#    with open("./flow.f90") as sourcefile:
#        sourcecode = sourcefile.read()
#        f2py.compile(sourcecode, modulename='flow', extension='.f90', extra_args=f_args)
#except:
#    print('Error in compileing flow!')
#    exit()



import glob
so_files = glob.glob('*.so')
print(so_files)
os.chdir('../')
from setuptools import setup
setup(name='util',
      version='0.1',
      description='Utility module',
      author='Yuchen Liu',
      author_email='yl8bc@mst.edu',
      license='MST',
      packages=['util'],
      data_files=[('util/', [f90_dir+filename for filename in so_files ])],
#                             './util/interp2d.so', \
#                             './util/uns.so', \
#                             './util/gridgen.so', \
#                             './util/mm.so', \
#                             './util/flow.so'])],
      zip_safe=False)
