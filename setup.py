import os.path
import setuptools
from glob import glob

requires = ["colorama",
            "pandas",
            "numpy",
            "matplotlib",
            "seaborn",
            "scikit-learn"]


with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
  name="saaps",
  version="1.3.1",
  author="Yang Xiao",
  author_email="fredrik1999@163.com",
  description="Single Amino Acids Polymorphism Statistics",
  keywords="SAP, Information Entropy, OneHot Encoding, Dimension Reduction",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/xiaosheep01/SAAPS",
  packages=setuptools.find_packages(),
  include_package_data=True,
  install_requires=requires,

  classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "License :: OSI Approved :: GNU Lesser General Public License v2 or later (LGPLv2+)",
        "Operating System :: OS Independent",
    ],
  entry_points={
             'console_scripts': [
                 'saaps = saaps.main:starts',
             ],
    }
)
