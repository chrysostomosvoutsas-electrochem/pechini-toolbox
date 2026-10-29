"""
PechiniLab — reproducible formulation calculator for Pechini-derived oxides.

Academic prototype, not a validated synthesis protocol.
The application calculates theoretical stoichiometry and reagent quantities from
explicit user-supplied precursor definitions. It does not predict phase purity,
reaction yield, pH, thermal behavior, or safe processing conditions.
"""

from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import re
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from typing import Dict, List, Mapping, Tuple

import pandas as pd
import streamlit as st

APP_NAME = "PechiniLab"
APP_VERSION = "0.2.0-academic-prototype"

# Standard atomic weights (g/mol), simplified fixed values for transparent calculations.
# For publication, document the source/version and consider standard uncertainty ranges.
ATOMIC_WEIGHTS: Dict[str, float] = {
    "H": 1.008, "He": 4.002602, "Li": 6.94, "Be": 9.0121831, "B": 10.81,
    "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998403163, "Ne": 20.1797,
    "Na": 22.98976928, "Mg": 24.305, "Al": 26.9815385, "Si": 28.085, "P": 30.973761998,
    "S": 32.06, "Cl": 35.45, "Ar": 39.948, "K": 39.0983, "Ca": 40.078,
    "Sc": 44.955908, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961, "Mn": 54.938044,
    "Fe": 55.845, "Co": 58.933194, "Ni": 58.6934, "Cu": 63.546, "Zn": 65.38,
    "Ga": 69.723, "Ge": 72.630, "As": 74.921595, "Se": 78.971, "Br": 79.904,
    "Kr": 83.798, "Rb": 85.4678, "Sr": 87.62, "Y": 88.90584, "Zr": 91.224,
    "Nb": 92.90637, "Mo": 95.95, "Ru": 101.07, "Rh": 102.9055, "Pd": 106.42,
    "Ag": 107.8682, "Cd": 112.414, "In": 114.818, "Sn": 118.710, "Sb": 121.760,
    "Te": 127.60, "I": 126.90447, "Xe": 131.293, "Cs": 132.90545196, "Ba": 137.327,
    "La": 138.90547, "Ce": 140.116, "Pr": 140.90766, "Nd": 144.242, "Sm": 150.36,
    "Eu": 151.964, "Gd": 157.25, "Tb": 158.92535, "Dy": 162.500, "Ho": 164.93033,
    "Er": 167.259, "Tm": 168.93422, "Yb": 173.045, "Lu": 174.9668, "Hf": 178.49,
    "Ta": 180.94788, "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217,
    "Pt": 195.084, "Au": 196.966569, "Hg": 200.592, "Tl": 204.38, "Pb": 207.2,
    "Bi": 208.98040, "Th": 232.0377, "U": 238.02891,
}

# Only simple element symbols, decimal coefficients, and nested () or [] groups are supported.
# Hydrate dot notation, charges, isotope notation, and phase labels are intentionally unsupported.
TOKEN_RE = re.compile(r"([A-Z][a-z]?|\d+(?:\.\d+)?|\.\d+|[()\[\]])")


class FormulaError(ValueError):
    """Raised when a chemical formula is invalid or unsupported."""


def parse_formula(formula: str) -> Dict[str, float]:
    """Parse a neutral inorganic formula with nested parentheses/brackets and decimals."""
    if not isinstance(formula, str) or not formula.strip():
        raise FormulaError("Enter a non-empty chemical formula.")
    source = re.sub(r"\s+", "", formula)
    tokens = TOKEN_RE.findall(source)
    if "".join(tokens) != source:
        raise FormulaError(
            "Unsupported formula syntax. Use element symbols, numeric coefficients, "
            "and balanced () or [] groups. Hydrate dots and charges are not supported."
        )

    index = 0

    def parse_group(expected_close: str | None = None) -> Dict[str, float]:
        nonlocal index
        result: Dict[str, float] = {}
        while index < len(tokens):
            token = tokens[index]
            if token in (")", "]"):
                if expected_close is None or token != expected_close:
                    raise FormulaError(f"Unexpected closing bracket '{token}'.")
                index += 1
                return result
            if token in ("(", "["):
                opening = token
                closing = ")" if opening == "(" else "]"
                index += 1
                nested = parse_group(closing)
                multiplier = 1.0
                if index < len(tokens) and re.fullmatch(r"(?:\d+(?:\.\d+)?|\.\d+)", tokens[index]):
                    multiplier = float(tokens[index])
                    index += 1
                if multiplier <= 0:
                    raise FormulaError("Group multipliers must be greater than zero.")
                for element, amount in nested.items():
                    result[element] = result.get(element, 0.0) + amount * multiplier
                continue
            if re.fullmatch(r"[A-Z][a-z]?", token):
                element = token
                if element not in ATOMIC_WEIGHTS:
                    raise FormulaError(f"Unknown element symbol '{element}'.")
                index += 1
                amount = 1.0
                if index < len(tokens) and re.fullmatch(r"(?:\d+(?:\.\d+)?|\.\d+)", tokens[index]):
                    amount = float(tokens[index])
                    index += 1
                if amount <= 0:
                    raise FormulaError("Element coefficients must be greater than zero.")
                result[element] = result.get(element, 0.0) + amount
                continue
            raise FormulaError(f"Unexpected token '{token}'.")
        if expected_close:
            raise FormulaError(f"Missing closing bracket '{expected_close}'.")
        return result

    composition = parse_group()
    if index != len(tokens):
        raise FormulaError("Formula was not fully parsed.")
    if not composition:
        raise FormulaError("No elements were found in the formula.")
    return composition


def molecular_weight(composition: Mapping[str, float]) -> float:
    """Return formula mass in g/mol using the fixed atomic weights above."""
    if not composition:
        raise FormulaError("Composition is empty.")
    unknown = set(composition) - set(ATOMIC_WEIGHTS)
    if unknown:
        raise FormulaError(f"Unknown element(s): {', '.join(sorted(unknown))}")
    mass = sum(ATOMIC_WEIGHTS[el] * amount for el, amount in composition.items())
    if not math.isfinite(mass) or mass <= 0:
        raise FormulaError("Calculated formula mass is not physically meaningful.")
    return mass


@dataclass
class Precursor:
    """Explicit precursor definition; MW is for the as-supplied chemical formula."""
    element: str
    name: str
    formula: str
    mw_g_mol: float
    target_atoms_per_formula: float
    assay_wt_percent: float
    phase: str = "solid"
    density_g_ml: float | None = None
    assay_basis: str = "active compound mass fraction"
    reference: str = "User supplied; not independently verified"

    def validate(self) -> None:
        if self.element not in ATOMIC_WEIGHTS:
            raise ValueError(f"Unknown target element: {self.element}")
        if not self.name.strip():
            raise ValueError("Precursor name is required.")
        if not self.formula.strip():
            raise ValueError(f"Formula is required for {self.name}.")
        parsed = parse_formula(self.formula)
        if self.element not in parsed:
            raise ValueError(
                f"Precursor formula {self.formula} does not contain target element {self.element}."
            )
        if self.target_atoms_per_formula <= 0:
            raise ValueError("Target-element atoms per precursor formula must be positive.")
        if not math.isfinite(self.mw_g_mol) or self.mw_g_mol <= 0:
            raise ValueError("Precursor molar mass must be positive.")
        if not (0 < self.assay_wt_percent <= 100):
            raise ValueError("Assay must be greater than 0 and at most 100 wt%.")
        if self.phase == "liquid" and (self.density_g_ml is None or self.density_g_ml <= 0):
            raise ValueError("A positive density is required for a liquid precursor.")
        if self.assay_basis != "active compound mass fraction":
            raise ValueError(
                "This version only supports assay reported as active-compound mass fraction."
            )


def calculate_precursor_amount(
    target_element_moles: float, precursor: Precursor
) -> Tuple[float, float, float]:
    """Return precursor moles, as-supplied mass (g), and optional liquid volume (mL)."""
    precursor.validate()
    if target_element_moles < 0 or not math.isfinite(target_element_moles):
        raise ValueError("Target element moles must be finite and non-negative.")
    precursor_moles = target_element_moles / precursor.target_atoms_per_formula
    active_mass = precursor_moles * precursor.mw_g_mol
    supplied_mass = active_mass / (precursor.assay_wt_percent / 100.0)
    volume_ml = (
        supplied_mass / precursor.density_g_ml if precursor.phase == "liquid" else float("nan")
    )
    return precursor_moles, supplied_mass, volume_ml


def calculate_batch(
    formula: str,
    product_mass_g: float,
    precursors: Mapping[str, Precursor],
    lithium_excess_percent: float = 0.0,
) -> dict:
    """Calculate theoretical elemental and precursor quantities; no yield/phase prediction."""
    composition = parse_formula(formula)
    mw = molecular_weight(composition)
    if not math.isfinite(product_mass_g) or product_mass_g <= 0:
        raise ValueError("Target theoretical product mass must be positive.")
    if not (0 <= lithium_excess_percent <= 100):
        raise ValueError("Optional lithium excess must be between 0 and 100%.")
    target_formula_moles = product_mass_g / mw
    element_rows = []
    precursor_rows = []
    for element, coefficient in composition.items():
        if element == "O":
            continue
        excess = 1.0 + lithium_excess_percent / 100.0 if element == "Li" else 1.0
        element_moles = target_formula_moles * coefficient * excess
        element_rows.append({
            "Element": element,
            "Atoms per formula unit": coefficient,
            "Target element amount (mol)": element_moles,
            "Applied excess (%)": lithium_excess_percent if element == "Li" else 0.0,
        })
        if element not in precursors:
            raise ValueError(
                f"No explicit precursor definition supplied for {element}. "
                "Add a precursor with a known formula, molar mass, assay basis, and physical properties."
            )
        precursor = precursors[element]
        pmol, mass_g, volume_ml = calculate_precursor_amount(element_moles, precursor)
        precursor_rows.append({
            "Element": element,
            "Precursor": precursor.name,
            "Precursor formula": precursor.formula,
            "Precursor amount (mol)": pmol,
            "Molar mass (g/mol)": precursor.mw_g_mol,
            "Assay (wt%)": precursor.assay_wt_percent,
            "Theoretical supplied mass (g)": mass_g,
            "Calculated volume (mL)": volume_ml,
            "Reference / provenance": precursor.reference,
        })
    return {
        "formula": formula,
        "product_mass_g": product_mass_g,
        "formula_mass_g_mol": mw,
        "target_formula_moles": target_formula_moles,
        "elements": element_rows,
        "precursors": precursor_rows,
        "lithium_excess_percent": lithium_excess_percent,
    }


def default_precursor_form(element: str) -> str:
    """Suggest a blank editable template only; does not infer a chemically valid precursor."""
    return ""


def main() -> None:
    st.set_page_config(page_title="PechiniLab | Academic Prototype", page_icon="⚗️", layout="wide")
    st.title("⚗️ PechiniLab")
    st.caption(
        "Academic prototype for transparent stoichiometry and precursor-quantity calculations. "
        "Not a validated synthesis protocol or safety assessment."
    )

    with st.sidebar:
        st.header("Calculation settings")
        st.write(f"Version: `{APP_VERSION}`")
        st.write("Atomic weights: fixed values embedded in this version.")
        st.warning(
            "Use only after checking the formula, precursor identity, assay basis, and "
            "independent calculations. No phase-purity, yield, pH, or thermal predictions are made."
        )
        st.markdown("**Supported formula syntax**")
        st.code("La0.6Sr0.4Co0.2Fe0.8O3\nLa(NO3)3\nLi7La3Zr2O12")
        st.caption("Hydrate dots, charges, isotope notation, and phase labels are not supported in target formulas.")

    col1, col2 = st.columns([1, 1])
    with col1:
        st.subheader("1. Target composition")
        formula = st.text_input("Target formula", value="LaCoO3", help="Use elemental coefficients; grouped formulas are supported.")
        product_mass = st.number_input(
            "Theoretical target product mass (g)", min_value=0.001, max_value=10000.0,
            value=10.0, step=0.5, format="%.3f"
        )
        li_excess = st.number_input(
            "Optional Li excess (wt% relative to stoichiometric Li amount)",
            min_value=0.0, max_value=100.0, value=0.0, step=1.0,
            help="This is a user-entered planning assumption, not a temperature-based prediction."
        )
        st.info(
            "Calculation basis: complete conversion to the specified nominal formula, "
            "with no allowance for yield, oxygen non-stoichiometry, phase impurities, or process losses."
        )

    try:
        composition = parse_formula(formula)
        formula_mass = molecular_weight(composition)
    except Exception as exc:
        st.error(f"Formula validation failed: {exc}")
        st.stop()

    with col2:
        st.subheader("2. Formula checks")
        st.metric("Formula mass", f"{formula_mass:.5f} g/mol")
        st.metric("Target formula amount", f"{product_mass / formula_mass:.6f} mol")
        st.dataframe(
            pd.DataFrame(
                [{"Element": el, "Coefficient": amount} for el, amount in composition.items()]
            ),
            use_container_width=True,
            hide_index=True,
        )

    st.divider()
    st.subheader("3. Define a precursor for each non-oxygen element")
    st.write(
        "All precursor data are explicit inputs. The app does not silently choose a default "
        "chemical form. Verify molecular masses and assay information against the supplier certificate "
        "or a citable reference before use."
    )

    precursors: Dict[str, Precursor] = {}
    element_cols = st.columns(min(max(len([e for e in composition if e != "O"]), 1), 3))
    for idx, element in enumerate([e for e in composition if e != "O"]):
        with element_cols[idx % len(element_cols)]:
            with st.container(border=True):
                st.markdown(f"**{element} precursor**")
                name = st.text_input(f"{element}: precursor name", key=f"name_{element}", placeholder="e.g., lanthanum nitrate hexahydrate")
                pformula = st.text_input(f"{element}: chemical formula", key=f"formula_{element}", placeholder="e.g., La(NO3)3")
                mw = st.number_input(
                    f"{element}: molar mass (g/mol)", min_value=0.001, value=1.0,
                    step=0.1, format="%.5f", key=f"mw_{element}"
                )
                atoms = st.number_input(
                    f"{element}: target atoms per precursor formula", min_value=0.001,
                    value=1.0, step=1.0, format="%.3f", key=f"atoms_{element}"
                )
                assay = st.number_input(
                    f"{element}: assay (wt%)", min_value=0.001, max_value=100.0,
                    value=99.0, step=0.5, format="%.3f", key=f"assay_{element}"
                )
                phase = st.selectbox(f"{element}: physical form", ["solid", "liquid"], key=f"phase_{element}")
                density = None
                if phase == "liquid":
                    density = st.number_input(
                        f"{element}: density (g/mL)", min_value=0.001, value=1.0,
                        step=0.01, format="%.4f", key=f"density_{element}"
                    )
                reference = st.text_input(
                    f"{element}: data source", key=f"ref_{element}",
                    placeholder="Supplier CoA, catalog number, DOI, or database version"
                )
                precursors[element] = Precursor(
                    element=element, name=name or f"UNSPECIFIED {element} precursor",
                    formula=pformula, mw_g_mol=float(mw),
                    target_atoms_per_formula=float(atoms), assay_wt_percent=float(assay),
                    phase=phase, density_g_ml=float(density) if density is not None else None,
                    reference=reference or "Not supplied; verify before use",
                )

    st.divider()
    st.subheader("4. Calculate and export")
    st.caption("The generated report is a calculation record, not an experimental SOP.")
    if st.button("Calculate theoretical reagent quantities", type="primary"):
        try:
            result = calculate_batch(formula, float(product_mass), precursors, float(li_excess))
        except Exception as exc:
            st.error(f"Calculation stopped: {exc}")
            st.stop()

        st.success("Calculation completed. Review every precursor definition and data source before relying on the values.")
        m1, m2, m3 = st.columns(3)
        m1.metric("Formula mass", f"{result['formula_mass_g_mol']:.5f} g/mol")
        m2.metric("Target formula amount", f"{result['target_formula_moles']:.6f} mol")
        m3.metric("Elements calculated", str(len(result["elements"])))

        st.markdown("#### Elemental stoichiometry")
        st.dataframe(pd.DataFrame(result["elements"]), use_container_width=True, hide_index=True)
        precursor_df = pd.DataFrame(result["precursors"])
        st.markdown("#### Precursor quantities")
        st.dataframe(
            precursor_df.style.format({
                "Precursor amount (mol)": "{:.6f}",
                "Molar mass (g/mol)": "{:.5f}",
                "Assay (wt%)": "{:.3f}",
                "Theoretical supplied mass (g)": "{:.5f}",
                "Calculated volume (mL)": lambda x: "—" if pd.isna(x) else f"{x:.4f}",
            }),
            use_container_width=True,
            hide_index=True,
        )

        record = {
            "application": APP_NAME,
            "application_version": APP_VERSION,
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "input": {
                "target_formula": formula,
                "theoretical_product_mass_g": float(product_mass),
                "optional_lithium_excess_percent": float(li_excess),
            },
            "calculation": result,
            "precursor_definitions": [asdict(p) for p in precursors.values()],
            "limitations": [
                "No experimental validation is implied.",
                "No phase purity, yield, pH, thermal behavior, or safety prediction is performed.",
                "Atomic weights are fixed constants embedded in this version.",
                "All precursor identities, molar masses, assay bases, and sources require independent verification.",
            ],
        }
        canonical = json.dumps(record, sort_keys=True, default=str, separators=(",", ":"))
        record["record_sha256"] = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
        st.download_button(
            "Download calculation record (JSON)",
            data=json.dumps(record, indent=2, default=str),
            file_name="pechinilab_calculation_record.json",
            mime="application/json",
        )
        st.download_button(
            "Download precursor table (CSV)",
            data=precursor_df.to_csv(index=False),
            file_name="pechinilab_precursor_quantities.csv",
            mime="text/csv",
        )
        st.markdown("#### Reproducibility and limitations")
        st.write(f"Record timestamp (UTC): {record['generated_utc']}")
        st.code(record["record_sha256"], language="text")
        st.warning(
            "Do not interpret these calculations as a recommendation to carry out a synthesis. "
            "Confirm the precursor chemistry, compatibility, hazards, and processing conditions "
            "with a qualified materials chemist and your institution's approved procedures."
        )

    with st.expander("Scientific scope and limitations"):
        st.markdown(
            """
            **Implemented:** formula parsing for a defined subset of formula syntax; formula mass;
            theoretical target formula moles; elemental amounts; precursor amounts based on an
            explicitly supplied number of target-element atoms per precursor formula; assay correction
            where assay means active-compound mass fraction; liquid volume from supplied density.

            **Not implemented:** chemical charge balance, oxidation-state inference, hydrate/solvate
            normalization, solution concentration models, thermodynamic speciation, complexation
            equilibria, pH prediction, polymerization kinetics, thermal analysis, lithium volatilization
            prediction, phase stability, yield prediction, uncertainty propagation, or experimental validation.

            **Publication status:** this is a software prototype. A research publication should report
            independent benchmark cases, software tests, data provenance, supported input syntax,
            limitations, versioned dependencies, and a clear statement of the scientific contribution.
            """
        )


if __name__ == "__main__":
    main()
