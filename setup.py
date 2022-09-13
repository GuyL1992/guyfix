from setuptools import setup, Extension


setup(name='spkmeans',
      version='1.0',
      description='implementation to calculate k clusters for given data points, for full or part the SPK process according to the user request',
      author= ' Guy Lamdan & Yair Ben Michael',
      ext_modules=[Extension('spkmeans', sources=['spkmeansmodule.c','spkmeans.c'])])
