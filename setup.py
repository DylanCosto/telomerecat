#!/usr/bin/env python3

from setuptools import Extension, setup
import os
import sys

extensions = [Extension("telomerecat._screening", ["telomerecat/_screening.c"], optional=True)]
# Opt in because this backend is compiled against the installed pysam headers.
# The default build keeps the portable public-API accelerator.
if os.environ.get('TELOMERECAT_BUILD_HTS_SCREENING') == '1':
  if not sys.platform.startswith('linux'):
    raise RuntimeError('The optional HTS screening build currently supports Linux')
  import pysam
  from Cython.Build import cythonize
  if pysam.config.HTSLIB != 'builtin':
    raise RuntimeError('The optional HTS screening build requires pysam with bundled HTSlib')
  library = next(path for path in pysam.get_libraries() if os.path.basename(path).startswith('libchtslib.'))
  extensions.extend(cythonize([Extension(
    'telomerecat._screening_hts', ['telomerecat/_screening_hts.pyx'],
    include_dirs=pysam.get_include(),
    define_macros=pysam.get_defines() + [('TELOMERECAT_PYSAM_VERSION', '"' + pysam.__version__ + '"')],
    library_dirs=[os.path.dirname(library)], libraries=[':' + os.path.basename(library)],
    runtime_library_dirs=['$ORIGIN/../pysam'],
  )], compiler_directives={'language_level': 3}))

setup(
  name="telomerecat",
  version='4.0.2',
  description="Telomere Computational Analysis Tool",
  long_description=open('README.md').read(),
  long_description_content_type='text/markdown',
  url='https://github.com/cancerit/telomerecat',
  author="JHR Farmery",
  license="GPL",
  python_requires='>= 3.7',
  author_email="cgphelp@sanger.ac.uk",
  packages=["telomerecat"],
  ext_modules=extensions,
  package_dir={"telomerecat": "telomerecat"},
  install_requires=["parabam>=3.0.1", "numpy", "pysam", "pandas", "click"],
  include_package_data=True,
  scripts=["./telomerecat/bin/telomerecat"],
  entry_points={'console_scripts': ['pysam_collate=telomerecat.pysam_collate:thin_wrap'],},
  zip_safe=False,
)
