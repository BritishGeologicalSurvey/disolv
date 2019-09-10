"""setup.py

Setup script for the disolv package

"""

from distutils.core import setup

setup(name='disolv',
      version='1.0',
      description='Python Distribution Utilities',
      author='Sarah Collins',
      author_email='sarcol@bgs.ac.uk',
      url='https://github.com/sarcol-bgs/disolv/',
      packages=["disolv"],
      install_requires=['numpy', 'scipy', 'matplotlib', 'pandas'],
      )
