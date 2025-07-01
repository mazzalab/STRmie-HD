from setuptools import setup, find_packages

setup(
    name="strmie",
    version="0.1.0",
    description="Tool per l'analisi di sequenze STR mitocondriali",
    author="Alessandro",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "pandas",
        "scikit-learn",
        "matplotlib",
        "tensorflow",
        "h5py",
        "joblib",
        "jinja2",
        "biopython",
        "seaborn"
    ],
    entry_points={
        "console_scripts": [
            "strmie=strmie.main:main"
        ]
    },
    include_package_data=True,
    zip_safe=False,
)
