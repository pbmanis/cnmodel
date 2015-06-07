from setuptools import setup
import os

def is_package(path):
    return (
        os.path.isdir(path) and
        os.path.isfile(os.path.join(path, '__init__.py'))
        )

def find_packages(path, base="" ):
    """ Find all packages in path """
    packages = {}
    for item in os.listdir(path):
        dir = os.path.join(path, item)
        if is_package( dir ):
            if base:
                module_name = "%(base)s.%(item)s" % vars()
            else:
                module_name = item
            packages[module_name] = dir
            packages.update(find_packages(dir, module_name))
    return packages
    
    
setup(name='nrnlibrary',
      version='0.2',
      description='Paul Manis Neuron Libraries in Python (Cells, Synapses, utilities)',
      url='',
      author='Paul Manis',
      author_email='pmanis@med.unc.edu',
      license='MIT',
      packages=find_packages("."),
      #data_files=[('nrnlibrary', ['x86_64/libnrnmech.la'])],
      zip_safe=False)
      
