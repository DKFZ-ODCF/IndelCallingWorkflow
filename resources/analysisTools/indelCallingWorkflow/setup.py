import ez_setup
ez_setup.use_setuptools()
from setuptools import setup

setup(
  name="biasFilter",
  version="0.1",
  author="Matthias Bieg",
  author_email="m.bieg@dkfz.de",
  description="Tool annotating PCR and sequencing biases to somatic mutations.",
  license="GPL verion 3",
  maintainer="Matthias Bieg",
  maintainer_email="m.bieg@dkfz.de",
  py_modules=["readbiasfunctions", "readbiasplots"],
  scripts=["biasFilter.py", "ez_setup.py"],
  install_requires=['matplotlib>=1.1.0', 'numpy>=1.6.2', 'pysam>=0.6', 'scipy>=0.11.0b1']
)
