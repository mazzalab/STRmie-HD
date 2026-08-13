from setuptools import setup, find_packages

setup(
    name="strmie",
    version="0.1.0",
    description="Tool for the analysis of CAG repeat sequences",
    author="Alessandro Napoli, Niccolò Liorni, Tommaso Mazza",
    packages=find_packages(),
    install_requires=['numpy', 'pandas', 'scikit-learn', 'matplotlib', 'h5py', 'joblib', 'jinja2', 'biopython', 'seaborn', 'colorama', 'openpyxl', 'xlrd', 'scipy'],
    extras_require={'bam': ['pysam']},
    entry_points={
        "console_scripts": [
            "strmie=strmie.main:main",
            "strmie-repeat=strmie.generic_main:main"
        ]
    },
    include_package_data=True,
    zip_safe=False,
)
