# STRmie

`STRmie` è un tool Python per l'analisi di dati biologici, con struttura a pacchetto e supporto per l'installazione moderna via `pyproject.toml`.

## ✅ Installazione (locale)

```bash
git clone https://github.com/mazzalab/STRmie.git
cd STRmie
conda env create -f searcHD.yml
conda activate searcHD
pip install -e .
```

## 🚀 Esecuzione

```bash
strmie
```

oppure:

```bash
python -m strmie.main
```

## 📁 Struttura

- `strmie/main.py`: file principale del tool
- `strmie/scripts/`: moduli e script personalizzati
- `model_object.pk`: modello salvato
- `searcHD.yml`: ambiente Conda (opzionale da aggiornare)
- `pyproject.toml`: configurazione moderna del pacchetto
