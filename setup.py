from setuptools import setup

setup(
    name="uptrop",
    version="2.0.0-alpha",
    packages = ['uptrop'],
    entry_points={
        'console_scripts': [
            'cloud_slice_tropomi=uptrop.main:main_tropomi',
        ],
    },
    author="Eloise A. Marais",
    author_email="e.marais@ucl.ac.uk.com",
    description="A Python package for processing TROPOMI satellite data",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/uptrop",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
    ],
    python_requires='>=3.7',
    install_requires=[
        "numpy",
        "matplotlib",
        "netCDF4",
    ],
)
