from setuptools import setup, Extension
import os
import platform

openmm_dir = '@OPENMM_DIR@'
exampleplugin_header_dir = '@EXAMPLEPLUGIN_HEADER_DIR@'
exampleplugin_library_dir = '@EXAMPLEPLUGIN_LIBRARY_DIR@'

# setup extra compile and link arguments on Mac
extra_compile_args = ['-std=c++11']
extra_link_args = []

if platform.system() == 'Darwin':
    extra_compile_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7']
    extra_link_args += ['-stdlib=libc++', '-mmacosx-version-min=10.7', '-Wl', '-rpath', openmm_dir+'/lib']

extension = Extension(name='_qtip4pforceplugin',
                      sources=['qTIP4PForcePluginWrapper.cpp'],
                      libraries=['OpenMM', 'qTIP4PForcePlugin'],
                      include_dirs=[os.path.join(openmm_dir, 'include'), exampleplugin_header_dir],
                      library_dirs=[os.path.join(openmm_dir, 'lib'), exampleplugin_library_dir],
                      extra_compile_args=extra_compile_args,
                      extra_link_args=extra_link_args
                     )

setup(name='qtip4pforceplugin',
      version='1.0',
      py_modules=['qtip4pforceplugin'],
      ext_modules=[extension],
     )
