# OpenGene Lab: Scientific Genetics Workbench & TUI

**OpenGene Lab** is a high-performance, modular Python CLI application and interactive terminal workbench designed for bioinformatics exploration, genetic drift simulation, and source-attributed genomic data management. Built for researchers, educators, and students, OpenGene Lab combines rigorous population genetics modeling with a professional terminal user interface.

---

##  Modular Architecture

The repository is structured into two core modules for clean separation of concerns:

* **`dataset.py`**: Manages data engineering, schema definitions, CSV generation across 12 genomic domains, provenance tracking, and text search querying.
* **`main.py`**: Houses the Wright-Fisher evolutionary simulation engine, cryptographic SHA-256 integrity verification, atomic file persistence, and the widescreen interactive TUI dashboard.

---

##  Key Features

### 1. Widescreen Scientific TUI (`btop`-Inspired)

* **High-Density Layout:** Engineered for widescreen terminals (1920x1200 viewports and adaptable standard terminals), featuring real-time Unicode sparkline frequency curves (` ▂▃▄▅▆▇█`), live telemetry, and multi-panel structured borders.
* **Interactive Menu System:** Seamlessly toggle between running new evolutionary simulations, querying local genomic datasets, verifying file integrity, and running system diagnostics.

### 2. Wright-Fisher Evolutionary Simulator

* **Neutral Genetic Drift Modeling:** Simulates independent neutral diploid Wright-Fisher replicates across generations with explicit reproduction assumptions.


* **Reproducible Seed Management:** Leverages NumPy's `SeedSequence.spawn` strategy to allocate isolated child random number generators per replicate, ensuring deterministic reproducibility.


* **Safety & Scalability Controls:** Built-in safeguards restrict excessive workloads (allele draws, population sizes, and generations) to prevent runaway resource consumption.



### 3. Source-Attributed Genomics Catalogue

* **Integrated Data Architecture:** Includes 12 standardized CSV datasets covering conditions, variants, structural abnormalities, phenotypes, genes, inheritance patterns, and population-level variation.
* **Explicit Provenance Tracking:** Preserves official database IDs, source URLs, retrieval dates, and conflicting clinical assertions (e.g., ClinVar submitter disagreements) instead of silently overwriting uncertain data.
* **Built-in CLI Query Engine:** Instantly search through cataloged genes, conditions, and phenotypes using case-insensitive substring matching across all columns.

### 4. Integrity & Security

* **Atomic File Operations:** Prevents file corruption during writes using temporary staging files and secure replacement mechanisms.


* **SHA-256 Checksum Verification:** Automatically generates and verifies cryptographic checksum manifests (`SHA256SUMS`) for all experiment run bundles to guarantee data integrity.



---

##  Quick Start Guide

### Prerequisites

* **Python 3.8+**
* **NumPy** (`python -m pip install numpy`)



### Installation

Clone the repository and ensure your working directory contains both `dataset.py` and `main.py`:

```bash
git clone https://github.com/x0ra-ssh/OpenGene-Lab.git
cd OpenGene-Lab
python -m pip install numpy

```

### Launching the Widescreen TUI

To open the interactive scientific terminal workbench, simply run:

```bash
python main.py

```

### Running a Simulation via CLI

To execute a neutral drift simulation with custom parameters and save a verified experiment bundle:

```bash
python main.py run --population-size 1000 --generations 200 --replicates 5 --seed 42 --output results/run-001

```

### Building & Querying Datasets

To generate the complete suite of source-attributed OpenGene Lab datasets (`opengene_data/`):

```bash
python main.py build-data

```

To query local datasets for a specific gene, phenotype, or condition:

```bash
python main.py query conditions --search "Marfan"

```

---

## ⚠️ Clinical & Academic Disclaimer

> **Disclaimer:** OpenGene Lab and its underlying simulation models and variant databases are strictly intended for **educational, academic research, and bioinformatics exploration purposes only**. They are **not** validated or approved for clinical diagnosis, patient prognosis, risk assessment, or individual medical treatment decisions.

---

## 📄 License

Distributed under the MIT License. See `LICENSE` for more information.
