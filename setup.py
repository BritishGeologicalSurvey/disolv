"""setup.py

Setup script for the disolv package

"""
from setuptools import setup

LONG_DESCRIPTION = """
# Disolv

> Modelling dilution tests.

See [https://github.com/sarcol-bgs/disolv/](https://github.com/sarcol-bgs/disolv/) for details.
"""

setup(name='disolv',
      version='1.4',
      license='GNU LGPLv3',
      description='Modelling dilution tests',
      long_description=LONG_DESCRIPTION,
      long_description_content_type='text/markdown',
      author='Sarah Collins',
      author_email='sarcol@bgs.ac.uk',
      url='https://github.com/sarcol-bgs/disolv/',
      download_url='https://github.com/sarcol-bgs/disolv/archive/v1.4.tar.gz',
      py_modules=['disolv', 'SolveEquation'],
      install_requires=['numpy', 'scipy', 'matplotlib', 'pandas'],
      classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 2.7",
        "License :: OSI Approved :: GNU General Public License v3 (LGPLv3)",
        "Operating System :: OS Independent",
      ],
      )
