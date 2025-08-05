#!/usr/bin/env python3
"""
REST2 Enhanced Sampling Project Setup
"""

from setuptools import setup, find_packages
import os

# Read README file
def read_readme():
    with open("README.md", "r", encoding="utf-8") as fh:
        return fh.read()

# Read requirements
def read_requirements():
    with open("requirements.txt", "r", encoding="utf-8") as fh:
        return [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="rest2-enhanced-sampling",
    version="1.0.0",
    author="REST2 Project Team",
    author_email="",
    description="Automated REST2 enhanced sampling simulation setup and execution",
    long_description=read_readme(),
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=read_requirements(),
    extras_require={
        "dev": [
            "pytest>=6.0.0",
            "pytest-cov>=3.0.0",
            "black>=22.0.0",
            "flake8>=4.0.0",
        ],
        "plotting": [
            "matplotlib>=3.5.0",
            "seaborn>=0.11.0",
        ],
        "jupyter": [
            "jupyter>=1.0.0",
            "ipykernel>=6.0.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "rest2-setup=main:main",
        ],
    },
    include_package_data=True,
    zip_safe=False,
) 