# ⚗️ PechiniLab: Precision Sol-Gel Engine for Electroceramics

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://pechini-toolbox.streamlit.app)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)

An open-source, physical-chemistry-grounded planning suite and automated bench protocol generator for the **Pechini (polymeric precursor) synthesis** of advanced electroceramics. 

Built specifically for materials in:
- **Solid-State Lithium Batteries (SSB):** Fast ion conductors ($\text{Li}_{0.4}\text{La}_{0.6}\text{TiO}_3$ [LLTO], $\text{Li}_7\text{La}_3\text{Zr}_2\text{O}_{12}$ [LLZO], $\text{Li}_{1+x}\text{Al}_x\text{Ti}_{2-x}(\text{PO}_4)_3$ [LATP]) and high-voltage spinels ($\text{LiNi}_{0.5}\text{Mn}_{1.5}\text{O}_4$ [LNMO]).
- **Solid Oxide Fuel Cells & Electrolyzers (SOFC / SOEC):** Mixed ionic-electronic conductors ($\text{La}_{1-x}\text{Sr}_x\text{Co}_{1-y}\text{Fe}_y\text{O}_{3-\delta}$ [LSCF], $\text{LaCo}_{0.8}\text{Fe}_{0.2}\text{O}_3$, $\text{La}_{1-x}\text{Sr}_x\text{MnO}_3$ [LSM]) and protonic electrolytes ($\text{BaZr}_{0.8}\text{Y}_{0.2}\text{O}_{3-\delta}$ [BZY]).
- **Electrocatalysts (OER / AEM / PEM):** Catalytic spinels and perovskites ($\text{NiFe}_2\text{O}_4$, $\text{Co}_3\text{O}_4$, $\text{LaNiO}_3$).

---

## 🔬 Scientific Rationale & Sol-Gel Physical Chemistry

Standard stoichiometric calculators treat sol-gel routes as simple solid-state powder mixtures ($m = n \times M$). In real-world polymeric precursor syntheses, phase purity and batch reproducibility fail due to **speciation equilibria, premature alkoxide hydrolysis, and foam overflow during charring**. 

**PechiniLab** implements governing physical and coordination chemistry principles to prevent these laboratory failure modes:

### 1. Coordination Thermodynamics & Valency-Driven Chelation ($R_{\text{CA/M}}$)
Citric acid ($\text{H}_4\text{Cit}$) coordinates metal cations ($M^{z+}$) via multidentate chelate rings. Low-valence cations ($\text{Li}^+$, $\text{Sr}^{2+}$, $\text{Ba}^{2+}$) exhibit weak stability constants ($\log K_1 \approx 0.5\text{–}3.0$), whereas transition and actinide/lanthanide cations ($\text{La}^{3+}$, $\text{Fe}^{3+}$, $\text{Ti}^{4+}$, $\text{Zr}^{4+}$) form strong chelates ($\log K_1 > 8\text{–}18$). 
- The engine dynamically calculates the minimum **Citric Acid to Total Metal Cation ratio ($R_{\text{CA/M}}$)** based on the average cation formal charge ($\bar{z}$):
  $$R_{\text{CA/M}} = \max\left(1.5,\, 1.0 + 0.5 \cdot \bar{z}\right)$$
- Transition metal alkoxides ($\text{Ti}(\text{OBu})_4$, $\text{Zr}(\text{OPr})_4$) trigger an automatic **hydrolysis-suppression regime** ($R_{\text{CA/M}} \ge 3.0\text{–}3.5$) to prevent localized condensation into hydrous oxide precipitates ($\text{TiO}_2 \cdot x\text{H}_2\text{O}$).

### 2. Speciation & Strict Precipitation-Free pH Windows
- Citrate carboxylate deprotonation ($\text{p}K_{\text{a}1} = 3.13$, $\text{p}K_{\text{a}2} = 4.76$, $\text{p}K_{\text{a}3} = 6.40$) requires base neutralization ($\text{NH}_4\text{OH}$) to complete metal coordination.
- The engine enforces distinct, domain-specific pH boundaries:
  - **SOFC / Nitrate Systems:** Strict target $\text{pH} = 5.5\text{–}6.5$. Exceeding $\text{pH} = 7.5$ triggers irreversible precipitation of divalent transition metal hydroxides ($\text{Co(OH)}_2$, $\text{Fe(OH)}_3$).
  - **Alkoxide Systems:** Target $\text{pH} = 7.0\text{–}7.5$ for alkoxide exchange without acid-catalyzed condensation.

### 3. Flory-Stockmayer Gelation Window ($R_{\text{EG/CA}}$)
- Transesterification between citric acid ($f=3$) and ethylene glycol ($f=2$) establishes a rigid, isotropic 3D polyester network. 
- The engine optimizes the crosslink density:
  - **Low-density / catalytic powders (SOFC):** $\text{EG:CA} \approx 2.0\text{–}2.5$ for clean, low-temperature combustion.
  - **High-density dense electrolytes (LLTO, LLZO):** $\text{EG:CA} \approx 4.0\text{–}8.0$ to mechanically entrap volatile ions.

### 4. Dynamic Precursor Assay & Loss on Ignition (LOI) Calibration
Hydrated nitrates ($\text{La(NO}_3)_3\cdot 6\text{H}_2\text{O}$, $\text{Fe(NO}_3)_3\cdot 9\text{H}_2\text{O}$) are deliquescent; moisture absorption routinely shifts effective assay by $2\text{–}5\%$. PechiniLab incorporates an **Assay (%) input** per reagent, recalculating raw weighed masses to maintain exact cation stoichiometry.

### 5. High-Temperature Volatilization Physics
Lithium oxide ($\text{Li}_2\text{O}$) undergoes sublimation via Knudsen-Langmuir evaporation during calcination above $750^\circ\text{C}$. The engine calculates temperature-dependent excess cation compensation (up to $+15\%$ for garnets at $1100^\circ\text{C}$) to guarantee single-phase product retention.

### 6. Crucible Sizing & Foam Expansion Safeguard
During the $200\text{–}250^\circ\text{C}$ primary charring stage, polyester resins expand $5\times\text{–}8\times$ in volume as gaseous $\text{H}_2\text{O}$, $\text{CO}_2$, and $\text{NO}_x$ evolve. PechiniLab calculates and enforces a **Minimum Crucible Volume ($\ge 35\ \text{mL/g}$ oxide)** to prevent furnace coil overflow.

---

## 🚀 Key Features

- **Arbitrary Formula Parsing:** Input any complex stoichiometry (e.g., `La0.6Sr0.4Co0.2Fe0.8O3`, `Li0.40La0.60TiO3`, `NiFe2O4`).
- **Comprehensive Precursor Database:** Out-of-the-box support for metal nitrates, acetates, chlorides, and liquid alkoxides ($\text{Ti}(\text{OBu})_4$, $\text{Zr}(\text{OPr})_4$).
- **Universal Cation Fallback Engine:** Dynamic estimation of atomic weights, standard oxidation states, and typical hydration states for unlisted or exotic dopants ($\text{Pr}$, $\text{Gd}$, $\text{Nb}$, $\text{Ta}$, $\text{Sm}$).
- **Step-by-Step Bench SOP:** Outputs a timestamped checklist with instructions for sequential alkoxide solvolysis, exotherm-controlled ammonia buffering, and multi-stage furnace heating profiles.

---

## 🛠️ Installation & Local Usage

To run the application locally:

```bash
# 1. Clone the repository
git clone [https://github.com/](https://github.com/)<your-username>/pechini-toolbox.git
cd pechini-toolbox

# 2. Install dependencies
pip install -r requirements.txt

# 3. Launch the Streamlit application
streamlit run app.py
