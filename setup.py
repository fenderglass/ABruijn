#(c) 2016-2019 by Authors
#This file is a part of the Flye package
#Released under the BSD license (see LICENSE file)

import os
import sys
import subprocess
import shutil

try:
    import setuptools
except ImportError:
    sys.exit("setuptools package not found. "
             "Please use 'pip install setuptools' first")

from setuptools import setup
from setuptools.command.install import install as SetuptoolsInstall
from distutils.command.build import build as DistutilsBuild
from distutils.spawn import find_executable

# Make sure we're running from the setup.py directory.
script_dir = os.path.dirname(os.path.realpath(__file__))
if script_dir != os.getcwd():
    os.chdir(script_dir)

from flye.__version__ import __version__


class MakeBuild(DistutilsBuild):
    def run(self):
        DistutilsBuild.run(self)

        if not find_executable("make"):
            sys.exit("ERROR: 'make' command is unavailable")
        try:
            subprocess.check_call(["make"])
        except subprocess.CalledProcessError as e:
            sys.exit("Compilation error: ", e)


class MakeInstall(SetuptoolsInstall):
    def run(self):
        SetuptoolsInstall.run(self)

        print('Installing C++ binaries')
        if os.path.isdir(self.install_lib) and not os.access(self.install_lib, os.W_OK):
            sys.exit('Error: no write permission for ' + self.install_lib + '  ' +
                     'Perhaps you need to use sudo?')

        if os.path.isdir(self.install_scripts) and not os.access(self.install_scripts, os.W_OK):
            sys.exit('Error: no write permission for ' + self.install_scripts + '  ' +
                     'Perhaps you need to use sudo?')

        if not os.path.exists(self.install_scripts):
            os.makedirs(self.install_scripts)

        build_dir = os.path.join(script_dir, "bin")
        install_dir = self.install_scripts
        bin_files = ['flye-modules', 'flye-minimap2', 'flye-samtools']
        for file in bin_files:
            if not os.path.isfile(os.path.join(build_dir, file)):
                sys.exit('Error: binary not found: ' + file)
            shutil.copy2(os.path.join(build_dir, file),
                         os.path.join(install_dir, file))


setup(name='flye',
      version=__version__,
      description='De novo assembler for single molecule sequencing reads using repeat graphs',
      url='https://github.com/fenderglass/Flye',
      author='Mikhail Kolmogorov',
      author_email = 'fenderglass@gmail.com',
      license='BSD-3-Clause',
      packages=['flye', 'flye/assembly', 'flye/config', 'flye/polishing',
                'flye/utils', 'flye/repeat_graph', 'flye/short_plasmids',
                'flye/trestle', 'flye/tests'],
      package_data={'flye': ['config/bin_cfg/*', 'tests/data/*']},
      entry_points={'console_scripts': ['flye = flye.main:main']},
      cmdclass={'build': MakeBuild,
                'install' : MakeInstall}
      )
