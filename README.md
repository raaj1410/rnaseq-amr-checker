# RNA-seq + AMR Conservative Checker (BV-BRC + NCBI SRA + ENA)

A local Python tool that checks, for one or more bacterial species:

1) **BV-BRC genome metadata** (genome_id, genome_name, strain, taxon_id, BioSample, BioProject)  
2) **BV-BRC AMR phenotype data** (filtered to lab-only conservative evidence)  
3) **RNA-seq availability** using **BioSample-only strong matching** from:
   - **NCBI SRA** (runinfo via E-utilities)
   - **ENA** (filereport via ENA portal API)

This tool is designed to minimise false positives. It may return fewer matches by design.

---
**Quick start** if you're in a hurry:
Windows (PowerShell)
```powershell
python -m venv .venv
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
pip install -r requirements.txt
python rnaseq_amr_checker.py --species "Pseudomonas aeruginosa"
```
macOS / Linux
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python rnaseq_amr_checker.py --species "Pseudomonas aeruginosa"
```
---

## What “conservative” means here
- AMR records are filtered to keep **lab-only** evidence and reject “predicted / in silico / computational” signals.
- RNA-seq is marked available only when a genome’s **BioSample accession** maps to real RNA-seq run IDs (**SRR/ERR/DRR**) using SRA and/or ENA.

---

## Input (how you provide species)

You can provide input in 3 ways:

### Option A: `--species` (comma-separated)
Single species:
```bash
python rnaseq_amr_checker.py --species "Pseudomonas aeruginosa"
```

Multiple species:
```bash
python rnaseq_amr_checker.py --species "Enterococcus faecium, Staphylococcus aureus"
```

### Option B: `--species-file` (one species per line)
Create a file `species.txt`:
```txt
Pseudomonas aeruginosa
Acinetobacter baumannii
Klebsiella pneumoniae
```

Run:
```bash
python rnaseq_amr_checker.py --species-file species.txt
```

### Option C: no input
If you run without `--species` or `--species-file`, the script uses the built-in default list inside the code.

---

## Output rules (important)

### Per-species outputs (always created)
For each species that succeeds:
- `SPECIES_LABONLY_conservative_rnaseq_amr.csv`
- `SPECIES_debug_stats.csv`
- `SPECIES_amr_taxon_summary.csv`

### Single-species run behaviour
If **exactly 1 species** succeeds:
- You get **CSV outputs only**
- **No ZIP**
- **No combined file**

### Multi-species run behaviour
If **2+ species** succeed:
- You get all per-species CSVs
- You also get:
  - `rnaseq_true_combined_<timestamp>.csv`
  - `rnaseq_amr_outputs_<timestamp>.zip` (contains all outputs including the combined CSV)

---

## RNASEQ availability column (you must filter it)
In the main output CSV (`*_LABONLY_conservative_rnaseq_amr.csv`), the column:

- `rnaseq_available` contains **both True and False**

If you want only RNA-seq-positive genomes, filter where:
- `rnaseq_available == True`

### Excel
Filter the `rnaseq_available` column → select `TRUE`

### Python (pandas)
```python
import pandas as pd
df = pd.read_csv("Pseudomonas_aeruginosa_LABONLY_conservative_rnaseq_amr.csv")
df_true = df[df["rnaseq_available"] == True]
```

---

## Requirements
- Python **3.9+** recommended (3.10 or 3.11 is ideal)
- Packages:
  - pandas
  - requests
  - tqdm
  - pyarrow

---

## Setup (recommended)

### Windows (PowerShell)

#### 1) Create venv
```powershell
python -m venv .venv
```

#### 2) If activation is blocked (common), run this and activate
Copy-paste both lines:
```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
```

#### 3) Install dependencies
If you have `requirements.txt`:
```powershell
pip install -r requirements.txt
```

Or install directly:
```powershell
pip install pandas requests tqdm pyarrow
```

---

### macOS/Linux

#### 1) Create + activate venv
```bash
python3 -m venv .venv
source .venv/bin/activate
```

#### 2) Install dependencies
```bash
pip install -r requirements.txt
```

Or:
```bash
pip install pandas requests tqdm pyarrow
```

---

## If venv creation fails (override options)

### Windows: Use Python Launcher
If `python` points to the wrong installation:
```powershell
py -3 -m venv .venv
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
```

### Windows: Rebuild a broken venv
```powershell
rmdir /s /q .venv
py -3 -m venv .venv
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
```

### Linux: Install missing venv support
If you see errors like `ensurepip is not available`:
```bash
sudo apt-get update
sudo apt-get install -y python3-venv python3-pip
python3 -m venv .venv
source .venv/bin/activate
```

### Last resort: Conda
```bash
conda create -n rnaseq-amr python=3.11 -y
conda activate rnaseq-amr
pip install pandas requests tqdm pyarrow
```

---

## NCBI environment variables (optional but recommended)
NCBI E-utilities rate limits are better with an API key.

### Windows (PowerShell)
```powershell
setx NCBI_EMAIL "your_email@example.com"
setx NCBI_API_KEY "your_ncbi_api_key"
```
Restart VS Code after setting these.

### macOS/Linux
```bash
export NCBI_EMAIL="your_email@example.com"
export NCBI_API_KEY="your_ncbi_api_key"
```

---

## Run the tool

### Default run (uses built-in species list)
```bash
python rnaseq_amr_checker.py
```

### Run a single species
```bash
python rnaseq_amr_checker.py --species "Pseudomonas aeruginosa"
```

### Run multiple species
```bash
python rnaseq_amr_checker.py --species "Acinetobacter baumannii, Klebsiella pneumoniae"
```

### Run from file
```bash
python rnaseq_amr_checker.py --species-file species.txt
```

### Custom output/cache folders
```bash
python rnaseq_amr_checker.py --cache-root "./cache" --out-root "./outputs"
```

---

## Output locations
- Outputs: `./outputs/`
- Cache: `./cache/`

Cache speeds up reruns. Delete `./cache/` if you want a clean re-fetch.

---

## Notes / limitations
- Input like `Enterobacter` (genus-only) may return 0 genomes because BV-BRC query expects an exact species string.
- Runtime depends heavily on:
  - number of BV-BRC taxon IDs processed
  - number of BioSamples after AMR filtering
  - API throttling / network stability
- Conservative filtering reduces matches by design.

---

## Disclaimer
This tool is for research/academic use. Always validate critical results with primary records and lab standards.
