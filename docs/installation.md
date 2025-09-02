# Installation

**Clone the repository:**
```bash
git clone https://github.com/mazzalab/STRmie-HD.git
cd STRmie-HD
```

**Create the conda environment:**
```bash
conda env create -f STRmie.yml
```

**Activate the environment:**
```bash
conda activate STRmie
```

**Install the package:**
```bash
pip install -e .
```

---

**Help**
```bash
strmie --help
```


**Automated Testing with Pytest**

The project includes a test suite that validates the core functionalities of both operational modes using example input and expected output files.  
This ensures the tool works as intended after installation or modification.

**Run tests with:**

```bash
pytest tests/test_strmie.py
```

**What it tests:**
Executes the full pipeline and compares the output Excel with the expected report.
