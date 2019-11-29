from setuptools import setup

setup(name='RNAlign2D',
      version='1.0.0',
      packages=['rnalign2d',],
      package_data={'': ['data/pseudo_matrix', 'data/simple_matrix']},
      include_package_data=True,
      url='',
      license='MIT',
      author='Tomasz Woźniak, Marcin Sajek',
      author_email='tomasz.wozniak@igcz.poznan.pl',
      description=(
          'Python tool allowing to use MUSCLE for secondary structure '
          'alignment of RNA'),
      entry_points={
          'console_scripts': [
              'rnalign2d = rnalign2d.rnalign2d:main',
              'create_matrix = rnalign2d.create_matrix:main',
              'rm_mod = rnalign2d.rm_mod:main'
          ]
      },
)
