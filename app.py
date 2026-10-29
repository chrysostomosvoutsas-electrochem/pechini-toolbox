"""
PechiniLab v0.3.0: Reproducible Synthesis Formulation Engine
Specialized for Perovskite (LLTO/LCTO) and Garnet (LLZO) Solid Electrolytes.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from typing import Dict, List, Mapping, Optional, Tuple

import pandas as pd
import streamlit as st

# Standard IUPAC atomic weights (g/mol)
ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.008, "Li": 6.94, "Be": 9.0122, "B": 10.81, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.9984, "Na": 22.9898, "Mg": 24.305, "Al": 26.9815,
    "Si": 28.085, "P": 30.9738, "S": 32.06, "Cl": 35.45, "K": 39.0983,
    "Ca": 40.078, "Sc": 44.9559, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961,
    "Mn": 54.9380, "Fe": 55.845, "Co": 58.9332, "Ni": 58.6934, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "Sr": 87.62, "Y": 88.9058,
    "Zr": 91.224, "Nb": 92.9064, "Mo": 95.95, "Ba": 137.327, "La": 138.9055,
    "Ce": 140.116, "Pr": 140.9077, "Nd": 144.242, "Sm": 150.36, "Eu": 151.964,
    "Gd": 157.25, "Ta": 180.9479, "W": 183.84, "Bi": 208.9804,
}

ORGANIC_REAGENTS = {
    "Citric Acid (anhydrous)": {"formula": "C6H8O7", "mw": 192.124, "density": 1.665},
    "Citric Acid (monohydrate)": {"formula": "C6H8O7*H2O", "mw": 210.14, "density": 1.54},
    "Ethylene Glycol": {"formula": "C2H6O2", "mw": 62.068, "density": 1.1132},
    "EDTA": {"formula": "C10H16N2O8", "mw": 292.24, "density": 0.86},
    "Urea": {"formula": "CH4N2O", "mw": 60.06, "density": 1.32},
}

PRECURSOR_LIBRARY: Dict[str, Dict[str, any]] = {
    "Li": {
        "name": "Lithium Nitrate",
        "formula": "LiNO3",
        "mw": 68.946,
        "atoms": 1,
        "phase": "solid",
        "handling": "Hygroscopic: dry at 120 °C prior to weighing to remove adsorbed water."
    },
    "La": {
        "name": "Lanthanum(III) Nitrate Hexahydrate",
        "formula": "La(NO3)3*6H2O",
        "mw": 433.01,
        "atoms": 1,
        "phase": "solid",
        "handling": "Confirm hydration by TGA or standardize. Alternatively use La2O3 calcined at 1000 °C dissolved in dilute HNO3."
    },
    "Ti": {
        "name": "Titanium(IV) Isopropoxide (TTIP)",
        "formula": "Ti(OCH(CH3)2)4",
        "mw": 284.22,
        "atoms": 1,
        "phase": "liquid",
        "density": 0.96,
        "handling": "Moisture-sensitive: pre-chelate in pure EtOH/dilute HNO3 with citric acid, or use peroxo-citrato route (Ti powder dissolved in H2O2 + NH3)."
    },
    "Zr": {
        "name": "Zirconyl Nitrate Hydrate",
        "formula": "ZrO(NO3)2*2H2O",
        "mw": 267.26,
        "atoms": 1,
        "phase": "solid",
        "handling": "For LLZO: verify ZrO2 residue wt% by calcining a test batch at 900 °C."
    },
    "Al": {
        "name": "Aluminum Nitrate Nonahydrate",
        "formula": "Al(NO3)3*9H2O",
        "mw": 375.13,
        "atoms": 1,
        "phase": "solid",
        "handling": "B-site dopant for LLTO/LLZO to stabilize cubic phase and enhance densification."
    },
    "Ta": {
        "name": "Tantalum(V) Chloride",
        "formula": "TaCl5",
        "mw": 358.21,
        "atoms": 1,
        "phase": "solid",
        "handling": "Reacts vigorously with water: dissolve in absolute ethanol containing citric acid."
    },
    "Sr": {
        "name": "Strontium Nitrate",
        "formula": "Sr(NO3)2",
        "mw": 211.63,
        "atoms": 1,
        "phase": "solid",
        "handling": "A-site dopant for LLTO to expand unit cell bottleneck and enhance conductivity."
    },
    "Eu": {
        "name": "Europium(III) Oxide",
        "formula": "Eu2O3",
        "mw": 351.926,
        "atoms": 2,
        "phase": "solid",
        "handling": "Dissolve in hot dilute HNO3 to generate europium nitrate in situ."
    }
}


def parse_sub_formula(sub_form: str) -> Dict[str, float]:
    """Parse single bracketed chemical formula group without hydrate dots."""
    token_re = re.compile(r"([A-Z][a-z]?|\d+(?:\.\d+)?|\.\d+|[()\[\]])")
    tokens = token_re.findall(sub_form.strip())
    if "".join(tokens) != sub_form.strip():
        raise ValueError(f"Invalid characters in sub-formula: {sub_form}")

    index = 0

    def parse_group(expected_close: Optional[str] = None) -> Dict[str, float]:
        nonlocal index
        result: Dict[str, float] = {}
        while index < len(tokens):
            token = tokens[index]
            if token in (")", "]"):
                if expected_close is None or token != expected_close:
                    raise ValueError(f"Unmatched closing bracket: {token}")
                index += 1
                return result
            if token in ("(", "["):
                closing = ")" if token == "(" else "]"
                index += 1
                nested = parse_group(closing)
                mult = 1.0
                if index < len(tokens) and re.fullmatch(r"(?:\d+(?:\.\d+)?|\.\d+)", tokens[index]):
                    mult = float(tokens[index])
                    index += 1
                for el, val in nested.items():
                    result[el] = result.get(el, 0.0) + val * mult
                continue
            if re.fullmatch(r"[A-Z][a-z]?", token):
                el = token
                if el not in ATOMIC_WEIGHTS:
                    raise ValueError(f"Unknown element: {el}")
                index += 1
                cnt = 1.0
                if index < len(tokens) and re.fullmatch(r"(?:\d+(?:\.\d+)?|\.\d+)", tokens[index]):
                    cnt = float(tokens[index])
                    index += 1
                result[el] = result.get(el, 0.0) + cnt
                continue
            raise ValueError(f"Syntax error near '{token}'")
        if expected_close:
            raise ValueError(f"Expected closing bracket: {expected_close}")
        return result

    return parse_group()


def parse_formula(formula: str) -> Dict[str, float]:
    """Parse chemical formulas supporting hydrate notation (e.g., La(NO3)3*6H2O or ·6H2O)."""
    formula = re.sub(r"\s+", "", formula)
    parts = re.split(r"[·*]", formula)
    total_comp: Dict[str, float] = {}

    for i, part in enumerate(parts):
        multiplier = 1.0
        m = re.match(r"^(\d+(?:\.\d+)?)(.*)$", part)
        if m and i > 0:  # e.g., 6H2O
            multiplier = float(m.group(1))
            part = m.group(2)
        sub_comp = parse_sub_formula(part)
        for el, cnt in sub_comp.items():
            total_comp[el] = total_comp.get(el, 0.0) + cnt * multiplier

    return total_comp


def molecular_weight(composition: Mapping[str, float]) -> float:
    return sum(ATOMIC_WEIGHTS[el] * amt for el, amt in composition.items())


@dataclass
class PrecursorConfig:
    element: str
    name: str
    formula: str
    mw_g_mol: float
    target_atoms_per_formula: float
    assay_wt_percent: float
    phase: str = "solid"
    density_g_ml: Optional[float] = None
    handling_notes: str = ""


def calculate_pechini_synthesis(
    nominal_formula: str,
    target_mass_g: float,
    precursor_map: Dict[str, PrecursorConfig],
    excess_map: Dict[str, float],
    ca_to_metal_ratio: float = 2.0,
    eg_to_ca_ratio: float = 4.0,
    urea_per_mass_oxide: float = 0.0,
) -> dict:
    comp = parse_formula(nominal_formula)
    mw_product = molecular_weight(comp)
    moles_product = target_mass_g / mw_product

    elements_calc = []
    precursors_calc = []
    total_metal_moles = 0.0

    for el, coeff in comp.items():
        if el == "O":
            continue
        applied_excess = excess_map.get(el, 0.0)
        target_moles = moles_product * coeff * (1.0 + applied_excess / 100.0)
        total_metal_moles += target_moles

        p = precursor_map[el]
        p_moles = target_moles / p.target_atoms_per_formula
        p_active_mass = p_moles * p.mw_g_mol
        p_supplied_mass = p_active_mass / (p.assay_wt_percent / 100.0)
        p_vol = (
            p_supplied_mass / p.density_g_ml
            if p.phase == "liquid" and p.density_g_ml
            else None
        )

        elements_calc.append({
            "Element": el,
            "Formula Stoich": coeff,
            "Target Moles (mol)": target_moles,
            "Excess Applied (%)": applied_excess,
        })

        precursors_calc.append({
            "Element": el,
            "Precursor": p.name,
            "Formula": p.formula,
            "Precursor Moles (mol)": p_moles,
            "Molar Mass (g/mol)": p.mw_g_mol,
            "Supplied Mass (g)": p_supplied_mass,
            "Volume (mL)": p_vol if p_vol is not None else float("nan"),
            "Special Handling": p.handling_notes,
        })

    # Pechini Organic Phase Calculation
    ca_moles = total_metal_moles * ca_to_metal_ratio
    ca_mw = ORGANIC_REAGENTS["Citric Acid (anhydrous)"]["mw"]
    ca_mass_g = ca_moles * ca_mw

    eg_moles = ca_moles * eg_to_ca_ratio
    eg_mw = ORGANIC_REAGENTS["Ethylene Glycol"]["mw"]
    eg_density = ORGANIC_REAGENTS["Ethylene Glycol"]["density"]
    eg_mass_g = eg_moles * eg_mw
    eg_volume_ml = eg_mass_g / eg_density

    urea_mass_g = target_mass_g * urea_per_mass_oxide

    return {
        "nominal_formula": nominal_formula,
        "target_mass_g": target_mass_g,
        "mw_product": mw_product,
        "moles_product": moles_product,
        "total_metal_moles": total_metal_moles,
        "elements": elements_calc,
        "precursors": precursors_calc,
        "organics": {
            "ca_moles": ca_moles,
            "ca_mass_g": ca_mass_g,
            "eg_moles": eg_moles,
            "eg_mass_g": eg_mass_g,
            "eg_volume_ml": eg_volume_ml,
            "urea_mass_g": urea_mass_g,
        },
    }


def main():
    st.set_page_config(page_title="PechiniLab v0.3.0 | Solid Electrolyte Synthesizer", page_icon="⚗️", layout="wide")
    st.title("⚗️ PechiniLab: Advanced Ceramic Formulation Suite")
    st.markdown(
        "**Reproducible Sol-Gel & Polymerizable Complex Calculator** for Perovskites "
        "($Li_{3x}La_{2/3-x}TiO_3$, $Li_{0.3}Ca_{0.35}TaO_3$) and Garnets ($Li_7La_3Zr_2O_{12}$)."
    )

    with st.sidebar:
        st.header("1. Material Benchmark Presets")
        preset = st.selectbox(
            "Select Literature Preset",
            [
                "Custom Composition",
                "LLTO: Li0.33La0.55TiO3 (Perovskite - High Ionic Conduct.)",
                "Al-LLTO: Li0.33La0.55Ti0.99Al0.01O3 (Huang 2022)",
                "LLZO: Li7La3Zr2O12 (Cubic Garnet Solid Electrolyte)",
                "Al-LLZO: Li6.25Al0.25La3Zr2O12 (Stabilized Cubic Garnet)",
                "LCTO: Li0.30Ca0.35TaO3 (Mariappan 2016)",
            ]
        )

        preset_defaults = {
            "LLTO: Li0.33La0.55TiO3 (Perovskite - High Ionic Conduct.)": ("Li0.33La0.55TiO3", {"Li": 10.0}),
            "Al-LLTO: Li0.33La0.55Ti0.99Al0.01O3 (Huang 2022)": ("Li0.33La0.55Ti0.99Al0.01O3", {"Li": 10.0}),
            "LLZO: Li7La3Zr2O12 (Cubic Garnet Solid Electrolyte)": ("Li7La3Zr2O12", {"Li": 15.0}),
            "Al-LLZO: Li6.25Al0.25La3Zr2O12 (Stabilized Cubic Garnet)": ("Li6.25Al0.25La3Zr2O12", {"Li": 15.0}),
            "LCTO: Li0.30Ca0.35TaO3 (Mariappan 2016)": ("Li0.30Ca0.35TaO3", {"Li": 5.0}),
        }

        default_formula = "Li0.33La0.55TiO3"
        default_excess = {"Li": 10.0}
        if preset in preset_defaults:
            default_formula, default_excess = preset_defaults[preset]

        target_formula = st.text_input("Nominal Ceramic Formula", value=default_formula)
        batch_mass = st.number_input("Target Sintered Powder Mass (g)", min_value=0.1, value=5.0, step=1.0)

        st.subheader("2. Organic Polyester Matrix")
        ca_ratio = st.slider("Citric Acid : Total Cation (CA:M)", min_value=1.0, max_value=5.0, value=2.0, step=0.5)
        eg_ratio = st.slider("Ethylene Glycol : Citric Acid (EG:CA)", min_value=1.0, max_value=6.0, value=4.0, step=0.5)
        urea_mass_mult = st.number_input("Urea auxiliary fuel (g per g oxide, 0 for pure Pechini)", value=0.0, step=0.2)

    try:
        comp = parse_formula(target_formula)
        mw = molecular_weight(comp)
    except Exception as e:
        st.error(f"Formula parsing error: {e}")
        st.stop()

    non_oxygen_elements = [el for el in comp.keys() if el != "O"]

    st.subheader("3. Precursor Identification & Volatilization Allowances")
    cols = st.columns(len(non_oxygen_elements))
    precursor_map: Dict[str, PrecursorConfig] = {}
    excess_map: Dict[str, float] = {}

    for i, el in enumerate(non_oxygen_elements):
        lib_data = PRECURSOR_LIBRARY.get(el, {
            "name": f"{el} Precursor", "formula": el, "mw": ATOMIC_WEIGHTS[el],
            "atoms": 1, "phase": "solid", "density": None, "handling": "Standard analytical grade"
        })
        with cols[i]:
            with st.container(border=True):
                st.markdown(f"### **{el}** (Atoms: {comp[el]})")
                name = st.text_input(f"Name ({el})", value=lib_data["name"], key=f"name_{el}")
                form = st.text_input(f"Formula ({el})", value=lib_data["formula"], key=f"form_{el}")
                try:
                    auto_mw = molecular_weight(parse_formula(form))
                except Exception:
                    auto_mw = lib_data["mw"]
                p_mw = st.number_input(f"MW g/mol ({el})", value=float(auto_mw), format="%.3f", key=f"mw_{el}")
                atoms = st.number_input(f"{el} per precursor", value=float(lib_data.get("atoms", 1)), key=f"at_{el}")
                assay = st.number_input(f"Assay wt% ({el})", min_value=10.0, max_value=100.0, value=99.0, key=f"as_{el}")
                phase = st.selectbox(f"Phase ({el})", ["solid", "liquid"], index=0 if lib_data["phase"] == "solid" else 1, key=f"ph_{el}")
                density = None
                if phase == "liquid":
                    density = st.number_input(f"Density g/mL ({el})", value=float(lib_data.get("density", 1.0)), key=f"d_{el}")
                
                def_exc = default_excess.get(el, 0.0)
                vol_exc = st.number_input(f"Thermal Excess % ({el})", min_value=0.0, max_value=50.0, value=float(def_exc), key=f"ex_{el}")
                excess_map[el] = vol_exc
                
                handling = st.text_area(f"Handling / Notes", value=lib_data.get("handling", ""), key=f"hd_{el}", height=68)

                precursor_map[el] = PrecursorConfig(
                    element=el, name=name, formula=form, mw_g_mol=p_mw,
                    target_atoms_per_formula=atoms, assay_wt_percent=assay,
                    phase=phase, density_g_ml=density, handling_notes=handling
                )

    st.divider()

    calc = calculate_pechini_synthesis(
        target_formula, batch_mass, precursor_map, excess_map,
        ca_to_metal_ratio=ca_ratio, eg_to_ca_ratio=eg_ratio, urea_per_mass_oxide=urea_mass_mult
    )

    colA, colB = st.columns([2, 1])
    with colA:
        st.subheader("4. Calculated Batch Sheet")
        precursor_df = pd.DataFrame(calc["precursors"])
        st.dataframe(
            precursor_df[["Element", "Precursor", "Formula", "Supplied Mass (g)", "Volume (mL)", "Special Handling"]].style.format({
                "Supplied Mass (g)": "{:.4f}",
                "Volume (mL)": lambda x: "—" if pd.isna(x) else f"{x:.3f}",
            }),
            use_container_width=True,
            hide_index=True
        )

    with colB:
        st.subheader("5. Polyester Network Quantities")
        org = calc["organics"]
        st.metric("Total Metal Cations (Σ Mⁿ⁺)", f"{calc['total_metal_moles']:.5f} mol")
        st.metric("Anhydrous Citric Acid", f"{org['ca_mass_g']:.3f} g ({org['ca_moles']:.4f} mol)")
        st.metric("Ethylene Glycol", f"{org['eg_volume_ml']:.2f} mL ({org['eg_mass_g']:.3f} g)")
        if org["urea_mass_g"] > 0:
            st.metric("Auxiliary Urea Fuel", f"{org['urea_mass_g']:.3f} g")

    with st.expander("Recommended Multi-Stage Thermal & Processing Protocol", expanded=True):
        st.markdown(
            f"""
            1. **Chelation & Complexation ($60–90\,^\circ\text{{C}}$):** Dissolve metal salts in dilute nitric acid or chelant solution. 
               *If using Ti or Zr alkoxides*, perform peroxo-complexation ($H_2O_2/NH_3$) or pre-chelate in alcohol to prevent uncontrolled precipitation.
            2. **Esterification / Gel Formation ($120–150\,^\circ\text{{C}}$):** Slowly evaporate excess solvent under vigorous stirring until a high-viscosity polymeric resin forms.
            3. **Charring / Carbonization ($350–500\,^\circ\text{{C}}$, 2–4 h):** Decompose the organic matrix in air to yield porous black carbonaceous precursor powder.
            4. **Calcination ($800–900\,^\circ\text{{C}}$, 2–5 h):** Burn off residual carbon. In LLTO/LLZO, this forms submicron, phase-pure crystalline crystallites (~100 nm).
            5. **Sintering ($1150–1250\,^\circ\text{{C}}$):** Consolidate pellets. Cover pellets with sacrificial mother powder of identical composition to suppress volatile $Li_2O$ loss.
            """
        )

    st.subheader("6. Audit Trail & Provenance")
    record = {
        "app": "PechiniLab v0.3.0",
        "utc_time": datetime.now(timezone.utc).isoformat(),
        "input": {"formula": target_formula, "target_mass_g": batch_mass, "ca_ratio": ca_ratio, "eg_ratio": eg_ratio},
        "results": calc
    }
    canonical = json.dumps(record, sort_keys=True, default=str)
    sha256_hash = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    record["sha256"] = sha256_hash

    dcol1, dcol2 = st.columns([1, 1])
    dcol1.download_button(
        "Download Lab Record (JSON)",
        data=json.dumps(record, indent=2, default=str),
        file_name=f"pechini_{target_formula}_batch.json",
        mime="application/json"
    )
    dcol2.download_button(
        "Download Reagents (CSV)",
        data=pd.DataFrame(calc["precursors"]).to_csv(index=False),
        file_name=f"pechini_{target_formula}_precursors.csv",
        mime="text/csv"
    )
    st.caption(f"Reproducibility Hash: `{sha256_hash}`")


if __name__ == "__main__":
    main()
