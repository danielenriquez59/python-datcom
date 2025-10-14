"""
Setup script for PyDATCOM package.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README
readme_file = Path(__file__).parent / "pydatcom" / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""

setup(
    name="pydatcom",
    version="0.1.0",
    author="Converted from USAF Digital DATCOM",
    description="Python implementation of USAF Digital DATCOM aerodynamic analysis tool",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/pydatcom",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering",
        "License :: Public Domain",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.19.0",
        "pyyaml>=5.3.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=2.10.0",
        ],
        "plots": [
            "matplotlib>=3.3.0",
        ],
        "analysis": [
            "pandas>=1.1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "pydatcom=pydatcom.cli:main",
        ],
    },
    package_data={
        "pydatcom": [
            "config/*.yaml",
            "data/tables/*.yaml",
            "data/tables/*.json",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)

