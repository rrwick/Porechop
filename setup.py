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
from distutils.command.build import build
from distutils.core import Command
import subprocess
import multiprocessing
import fnmatch
import importlib.util

# Install setuptools if not already present.
if not importlib.util.find_spec("setuptools"):
    import ez_setup
    ez_setup.use_setuptools()

from setuptools import setup
from setuptools.command.install import install

# Get the program version from another file.
exec(open('porechop/version.py').read())

with open('README.md', 'rb') as readme:
    LONG_DESCRIPTION = readme.read().decode()


class PorechopBuild(build):
    """
    The build process runs the Makefile to build the C++ functions into a shared library.
    """

    def run(self):
        build.run(self)  # Run original build code

        clean_cmd = ['make', 'clean']
        try:
            make_cmd = ['make', '-j', str(min(8, multiprocessing.cpu_count()))]
        except NotImplementedError:
            make_cmd = ['make']

        def clean_cpp():
            subprocess.call(clean_cmd)

        def compile_cpp():
            subprocess.call(make_cmd)

        self.execute(clean_cpp, [], 'Cleaning previous compilation: ' + ' '.join(clean_cmd))
        self.execute(compile_cpp, [], 'Compiling Porechop: ' + ' '.join(make_cmd))


class PorechopInstall(install):
    """
    The install process copies the C++ shared library to the install location.
    """

    def run(self):
        install.run(self)  # Run original install code
        shutil.copyfile(os.path.join('porechop', 'cpp_functions.so'),
                        os.path.join(self.install_lib, 'porechop', 'cpp_functions.so'))


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
      cmdclass={'build': PorechopBuild,
                'install': PorechopInstall,
                'clean': PorechopClean}
      )
