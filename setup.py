"""
Setup script for OpenSNP Project 10

Install in development mode:
    pip install -e .

Install with dev dependencies:
    pip install -e ".[dev]"
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read README for long description
readme_file = Path(__file__).parent / "README.md"
long_description = readme_file.read_text() if readme_file.exists() else ""

setup(
    name="opensnp-project10",
    version="1.0.0",
    description="Genetic Similarity Network & Cryptic Relatedness in OpenSNP",
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="kmcerr",
    license="MIT",
    python_requires=">=3.10",
    packages=find_packages(exclude=["tests", "tests.*"]),
    install_requires=[
        "pandas>=1.5",
        "numpy>=1.23",
        "matplotlib>=3.6",
        "networkx>=2.8",
        "tqdm>=4.65",
    ],
    extras_require={
        "dev": [
            "pytest>=7.4",
            "pytest-cov>=4.1",
            "mypy>=1.5",
            "black>=23.0",
            "flake8>=6.0",
            "isort>=5.12",
        ],
    },
    entry_points={
        "console_scripts": [
            "opensnp-pipeline=run_pipeline:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
