#pip install pybind11
#python setup.py bdist_wheel

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup
import sys
import os
import subprocess

__version__ = "1.0.1"

#print('Compile siteloc')
rootdir = os.path.abspath(os.path.dirname(__file__))
tool_path = os.path.join(rootdir, 'msproteomics', 'tools')
print(tool_path)
p = subprocess.Popen(["make"], stdout = subprocess.PIPE, cwd = tool_path)
p.wait()
#print('Done.')

extra_compile_args = None
extra_link_args = None

# separate compiler options for Windows
if sys.platform.startswith('win'):
    extra_compile_args = ['-openmp']
# Use OpenMP if directed or not on a Mac
elif os.environ.get('USEOPENMP') or not sys.platform.startswith('darwin'):
    extra_compile_args = ['-fopenmp']
    extra_link_args = [
        '-lgomp'
    ]

ext_modules = [
    Pybind11Extension(
        "_msproteomics",
        ["src/iq_pybind11.cpp"],
        define_macros=[("VERSION_INFO", __version__), ('EIGEN_DONT_PARALLELIZE', None), ('_OPENMP', None)],
        extra_compile_args = extra_compile_args,
        extra_link_args = extra_link_args,
        include_dirs=["src/"]
    ),
]


setup(
    name = 'msproteomics',
    version = __version__,    
    description = 'A Python package for mass spectrometry-based proteomics data processing',
    url = 'https://github.com/tvpham/msproteomics',
    author = 'Thang Pham',
    author_email = 't.pham@amsterdamumc.nl',
    license = 'Apache License 2.0',
    packages = ['msproteomics'],
    long_description = """# msproteomics sitereport: reporting DIA-MS phosphoproteomics experiments at site level with ease\n""",
    long_description_content_type = 'text/markdown',
    classifiers = [
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
		'License :: OSI Approved :: Apache Software License',
    ],
    install_requires = [    
    'numpy',
    'pandas',
    'scipy',
    'matplotlib'
    ], 
    entry_points = {
        'console_scripts': ['sitereport=msproteomics.sitereport:main',
                            'read_diann=msproteomics.read_diann:main'],
    },
    ext_modules = ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    include_package_data=True,
    package_data={"": ["siteloc*"]}
)
