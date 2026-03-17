"""Setup configuration for the Hybrid Capture Probe Designer."""

from setuptools import setup, find_packages

setup(
    name="probe-designer",
    version="1.0.0",
    description=(
        "Hybrid Capture Probe Designer - Design ultra-short capture probes "
        "for DNA targets (<40 bp) with optional LNA incorporation."
    ),
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Probe Designer Contributors",
    license="MIT",
    packages=find_packages(),
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "probe-designer=probe_designer.cli:main",
        ],
    },
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
)
