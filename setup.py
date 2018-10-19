#!/usr/bin/env python3
"""
Run 'python3 setup.py install' to install Porechop.
"""

# Make sure this is being run with Python 3.4 or later.
import sys
if sys.version_info.major != 3 or sys.version_info.minor < 4:
    print('Error: you must execute setup.py using Python 3.4 or later')
    sys.exit(1)

import os
import shutil
import subprocess
import multiprocessing
import fnmatch
import importlib.util

# Install setuptools if not already present.
if not importlib.util.find_spec("setuptools"):
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup, Extension, Command
from setuptools.command.install import install

# Get the program version from another file.
exec(open('porechop/version.py').read())

with open('README.md', 'rb') as readme:
    LONG_DESCRIPTION = readme.read().decode()




class PorechopClean(Command):
    """
    Custom clean command that really cleans up everything, except for:
      - the compiled *.so file needed when running the programs
      - setuptools-*.egg file needed when running this script
    """
    user_options = []

    def initialize_options(self):
        self.cwd = None

    def finalize_options(self):
        self.cwd = os.getcwd()

    def run(self):
        assert os.getcwd() == self.cwd, 'Must be in package root: %s' % self.cwd

        delete_directories = []
        for root, dir_names, filenames in os.walk(self.cwd):
            for dir_name in fnmatch.filter(dir_names, '*.egg-info'):
                delete_directories.append(os.path.join(root, dir_name))
            for dir_name in fnmatch.filter(dir_names, 'build'):
                delete_directories.append(os.path.join(root, dir_name))
            for dir_name in fnmatch.filter(dir_names, '__pycache__'):
                delete_directories.append(os.path.join(root, dir_name))
        for delete_directory in delete_directories:
            print('Deleting directory:', delete_directory)
            shutil.rmtree(delete_directory)

        delete_files = []
        for root, dir_names, filenames in os.walk(self.cwd):
            for filename in fnmatch.filter(filenames, 'setuptools*.zip'):
                delete_files.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.o'):
                delete_files.append(os.path.join(root, filename))
            for filename in fnmatch.filter(filenames, '*.pyc'):
                delete_files.append(os.path.join(root, filename))
        for delete_file in delete_files:
            print('Deleting file:', delete_file)
            os.remove(delete_file)

#TODO: enable these
DEBUGFLAGS   = ['-DSEQAN_ENABLE_DEBUG=1', '-g']
RELEASEFLAGS = ['-O3', '-D', 'NDEBUG']
if 'debug' in sys.argv:
    release_args = ['-DSEQAN_ENABLE_DEBUG=1', '-g']
    sys.argv.remove('debug')
else:
    release_args = ['-O3', '-D', 'NDEBUG']

extensions = []

seqan_dirs = []
for dirpath, dirnames, _ in os.walk(os.path.join('porechop', 'include', 'seqan')):
    seqan_dirs.extend(os.path.join(dirpath, x) for x in dirnames)

extensions.append(Extension(
    'porechop.porechopHelpers',
    sources=[os.path.join('porechop', 'src', x)
        for x in ('adapter_align.cpp', 'alignment.cpp')
    ],
    include_dirs=[os.path.join('porechop', 'include')] + seqan_dirs,
    extra_compile_args=['-std=c++14', '-fPIC', '-Wall', '-Wextra', '-pedantic'] + release_args
))


setup(name='porechop',
      version=__version__,
      description='Porechop',
      long_description=LONG_DESCRIPTION,
      url='http://github.com/rrwick/porechop',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPL',
      packages=['porechop'],
      entry_points={"console_scripts": ['porechop = porechop.porechop:main']},
      zip_safe=False,
      ext_modules=extensions,
      cmdclass={'clean': PorechopClean},
      )
