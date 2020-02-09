from setuptools import setup

setup(name='espyranto',
      version='0.0.1',
      description='Module for multi-well plate data',
      url='http://github.com/espyranto/espyranto',
      author='John Kitchin',
      author_email='jkitchin@andrew.cmu.edu',
      license='GPL',
      platforms=[],
      packages=['espyranto'],
      scripts=[],
      include_package_data=True,
      long_description='''A module to read plate data.''',
      install_requires=['ase'],)
