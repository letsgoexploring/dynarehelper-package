from setuptools import setup

setup(
  name = 'dynarehelper',
  packages = ['dynarehelper'],
  version = '0.0.1',
  description = 'A package for using Dynare++ and Octave to solve dynamic stochastic general equilibrium (DSGE) models',
  author = 'Brian C. Jenkins',
  author_email = 'bcjenkin@uci.edu',
  url = 'https://github.com/letsgoexploring/dynarehelper-package',
  # download_url = 'https://github.com/letsgoexploring/dynarehelper-package/tarball/dist/dynarehelper-0.0.1.tar.gz',
  keywords = ['economics','solution','dsge','macroconomics','rbc','new keynesian','business cycles','dynare','dynare++','octave'],
  classifiers = [],
  package_data={'dynarehelper': ['dynare_simul_.mex','dynare_simul.m']},
  # package_data={'dynarehelper': ['dynare_simul_.mex','dynare_simul.m']},
  # data_files=[('.',['dynare_simul_.mex','dynare_simul.m'])],
  include_package_data=True,
)