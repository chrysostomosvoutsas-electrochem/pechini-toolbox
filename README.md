# PechiniLab: Solid Electrolyte Synthesis Formulation Engine

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)

**PechiniLab** is an open-source, reproducible formulation calculator tailored for sol-gel and polymerizable complex (Pechini) synthesis of advanced ceramic solid electrolytes:
- **Perovskites:** $Li_{3x}La_{2/3-x}TiO_3$ (LLTO) and $Li_{0.3}Ca_{0.35}TaO_3$ (LCTO)
- **Garnets:** $Li_7La_3Zr_2O_{12}$ (LLZO) and substituted variants (Al, Ga, Ta)
- **NASICON:** $Li_{1+x}Al_xTi_{2-x}(PO_4)_3$ (LATP)

---

## Key Features

1. **Hydrate-Aware Formula Parser:** Fully parses complex stoichiometry including bracketed ligands and hydrate dot notations (e.g., `La(NO3)3*6H2O`, `Al(NO3)3·9H2O`, `Ti(OCH(CH3)2)4`).
2. **Polyester Network Engine:** Computes exact stoichiometric quantities for chelating agents (Citric Acid) and cross-linking agents (Ethylene Glycol) based on literature cation-to-chelant molar ratios ($[CA]:[\Sigma M^{n+}]$ and $[EG]:[CA]$).
3. **Thermal Volatilization Allowance:** Custom element-by-element thermal excess percentages (e.g., +10% Li for LLTO, +15% Li for LLZO) to compensate for high-temperature calcination and sintering loss.
4. **Literature Presets:** Instant configuration based on published synthesis procedures (Vijayakumar 2004, Mariappan 2016, Huang 2022).
5. **FAIR Data & Provenance:** Exports cryptographically hashed (SHA-256) JSON laboratory records and CSV reagent tables for experimental notebooks.

---

## Installation & Usage

### 1. Clone repository
```bash
git clone https://github.com/your-username/pechinilab.git
cd pechinilab
```

### 2. Create and activate virtual environment
```bash
python -m venv .venv

# On Linux/macOS:
source .venv/bin/activate

# On Windows:
.venv\Scripts\activate
```

### 3. Install dependencies
```bash
pip install -r requirements.txt
```

### 4. Run the Streamlit Application
```bash
streamlit run app.py
```

### 5. Run Unit Tests
```bash
pytest -v
```

---

## Scientific Scope & Guidelines

- **Precursor Drying:** Reagents such as $LiNO_3$ and $Li_2CO_3$ must be oven-dried at 120 °C prior to analytical weighing. $La_2O_3$ should be pre-calcined at 1000 °C to eliminate adsorbed moisture and $La(OH)_3$.
- **Titanium/Zirconium Chelation:** Alkoxides hydrolyze violently in water. Always pre-chelate in alcohol or utilize the peroxo-citrato route ($Ti$ metal powder dissolved in aqueous $H_2O_2 + NH_3$).
- **Sintering Atmosphere:** Sintering lithium conductors above 1100 °C requires covering pressed pellets with sacrificial mother powder of identical composition to suppress lithium volatilization.

---

## Citation

If you use PechiniLab in your research, please cite:
```bibtex
@software{pechinilab2026,
  author = {Voutsas, Chrysostomos},
  title = {PechiniLab: Reproducible Formulation Suite for Solid-State Oxide Electrolytes},
  year = {2026},
  url = {https://github.com/your-username/pechinilab}
}
```

## License
MIT License.
