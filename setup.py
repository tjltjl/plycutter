from setuptools import setup, find_packages

setup(name='plycutter',
      version='0.0',
      description='Slice 3D shell models into finger-jointed laser-cut parts',
      url='http://github.com/tjltjl/plycutter',
      author='Tuomas Lukka',
      author_email='tuomas@hipcode.fi',
      license='AGPL-3.0-or-later',
      python_requires='>=3.7.0',
      packages=find_packages(),
      zip_safe=True,
      tests_require=[
          'hypothesis',
          'pytest-benchmark',
      ],
      install_requires=[
          'numpy',
          'pyrsistent',
          'ezdxf',
          'gmpy2',
          'trimesh',
          'shapely',
          'matchpy',
          'sortedcontainers',
          'matplotlib',
          'mkdocs',
          'mkdocs-material',
          'mkdocstrings',
          'scipy',  # Soft dep of part of trimesh we need
          'networkx',  # Soft dep of part of trimesh we need
          'rtree',  # Soft dep of part of trimesh we need
      ],
      entry_points={
          'console_scripts': ['plycutter=plycutter.command_line:main'],
      },
      )
