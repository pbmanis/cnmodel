from setuptools import setup, find_packages
import os

path = os.path.join(os.path.dirname(__file__), 'cnmodel')
version = None
for line in open(os.path.join(path, '__init__.py'), 'r').readlines():
    if line.startswith('__version__'):
        version = line.partition('=')[2].strip('"\' \n')
        break
if version is None:
    raise Exception("Could not read __version__ from cnmodel/__init__.py")


setup(name='cnmodel',
      version='0.3',
      description='Neuron Libraries in Python (Cells, Synapses, anmodel, Decorator, Morphology, Populations, utilities)',
      url='',
      author='Paul B. Manis and Luke Campagnola',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages(include=['cnmodel*']),
      zip_safe=False)
      