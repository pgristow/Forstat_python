"""
Setup script for Forstat - Forensic Statistics Application
"""
from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = [line.strip() for line in fh if line.strip() and not line.startswith("#")]

setup(
    name="forstat",
    version="0.1.0",
    author="Forstat Development Team",
    description="Windows application for forensic and population genetics analysis",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/pgristow/Forstat_python",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: Microsoft :: Windows",
    ],
    python_requires=">=3.9",
    install_requires=requirements,
    entry_points={
        "console_scripts": [
            "forstat=main:main",
        ],
    },
    include_package_data=True,
    package_data={
        "forstat": [
            "resources/icons/*",
            "resources/themes/*",
            "resources/templates/*",
        ],
    },
)
