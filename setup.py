"""setup.py

Setup script for the disolv package

"""

from distutils.core import setup

setup(name='disolv',
      version='1.2',
	  licence='GNU LGPLv3'
      description='Modelling dilution tests',
      author='Sarah Collins',
      author_email='sarcol@bgs.ac.uk',
      url='https://github.com/sarcol-bgs/disolv/',
      packages=["disolv"],
      install_requires=['numpy', 'scipy', 'matplotlib', 'pandas'],
      )
