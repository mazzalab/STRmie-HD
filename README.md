# STRmie

`STRmie` is a Python tool for biological data analysis, organized as a package and supporting modern installation via `pyproject.toml`.

## ✅ Installation (local)

```bash
git clone https://github.com/mazzalab/STRmie.git
cd STRmie
conda env create -f searcHD.yml
conda activate searcHD
pip install -e .
```

## 🚀 Execution

```bash
strmie
```

or:

```bash
python -m strmie.main
```

## 📁 Structure

- `strmie/main.py`: main entry point of the tool  
- `strmie/scripts/`: custom modules and scripts  
- `model_object.pk`: pre-trained model file  
- `searcHD.yml`: Conda environment file (optional, may need updating)  
- `pyproject.toml`: modern package configuration file  
