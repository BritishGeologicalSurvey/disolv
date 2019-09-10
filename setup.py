"""setup.py

Setup script for the disolv package

"""

from setuptools import find_packages, setup

setup(name='disolv',
      version='1.2',
      license='GNU LGPLv3',
      description='Modelling dilution tests',
      author='Sarah Collins',
      author_email='sarcol@bgs.ac.uk',
      url='https://github.com/sarcol-bgs/disolv/',
      download_url='https://github.com/sarcol-bgs/disolv/archive/v1.2.tar.gz',
      py_modules=['disolv', 'SolveEquation'],
      install_requires=['numpy', 'scipy', 'matplotlib', 'pandas'],
      classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU LGPLv3 License",
        "Operating System :: OS Independent",
      ],
      )
