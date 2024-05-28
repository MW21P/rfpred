from setuptools import setup, find_packages

setup(
    name="rfpred",
    version="0.1",
    author="Milena Wiegand and Matthias Galka",
    author_email="milena.wiegand@epfl.ch and matthias.galka@epfl.ch",
    description="A NN driven Rf prediction tool",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/MW21P/rfpred.git",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={"rfpred": ["data/*.csv"]},
   
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.10',
    install_requires=[
        "pandas",
        "numpy",
        "rdkit"
    ],
    extras_require=dict(
        test=[
            "pytest",
            "pytest-cov",
        ])
)
