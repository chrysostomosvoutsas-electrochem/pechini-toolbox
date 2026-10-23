import streamlit as st
import pandas as pd
import numpy as np
import re

# ---------------------------------------------------------
# PAGE SETUP & SCIENTIFIC STYLING
# ---------------------------------------------------------
st.set_page_config(
    page_title="PechiniLab | Universal Electroceramics Sol-Gel Suite",
    page_icon="⚗️",
    layout="wide"
)

st.title("⚗️ PechiniLab: Universal Sol-Gel Suite for Electrochemical Materials")
st.markdown(
    "A comprehensive stoichiometry & polyesterification engine for **Solid Oxide Fuel Cells (SOFC/SOEC)**, "
    "**Solid-State Batteries (SSB)**, and **Electrocatalysts (AEM/PEM/OER)**."
)

# ---------------------------------------------------------
# COMPREHENSIVE PERIODIC & REAGENT DATABASE
# ---------------------------------------------------------
ATOMIC_WEIGHTS = {
    "H": 1.008, "Li": 6.941, "Be": 9.012, "B": 10.81, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.998, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085,
    "P": 30.974, "S": 32.065, "Cl": 35.453, "K": 39.098, "Ca": 40.078, "Sc": 44.956,
    "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933,
    "Ni": 58.693, "Cu": 63.546, "Zn": 65.380, "Ga": 69.723, "Ge": 72.630, "Rb": 85.468,
    "Sr": 87.620, "Y": 88.906, "Zr": 91.224, "Nb": 92.906, "Mo": 95.950, "Ru": 101.07,
    "Rh": 102.91, "Pd": 106.42, "Ag": 107.87, "Cd": 112.41, "In": 114.82, "Sn": 118.71,
    "Sb": 121.76, "Te": 127.60, "Cs": 132.91, "Ba": 137.33, "La": 138.91, "Ce": 140.12,
    "Pr": 140.91, "Nd": 144.24, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93,
    "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Re": 186.21, "Ir": 192.22, "Pt": 195.08,
    "Au": 196.97, "Pb": 207.2, "Bi": 208.98
}

STANDARD_VALENCIES = {
    "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Ag": 1,
    "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "Zn": 2, "Ni": 2, "Cu": 2, "Co": 2, "Fe": 3, "Mn": 3, "Pb": 2,
    "Al": 3, "Sc": 3, "Y": 3, "La": 3, "Ga": 3, "Cr": 3, "Ce": 3, "Pr": 3, "Nd": 3, "Sm": 3, "Eu": 3, "Gd": 3, "In": 3, "Bi": 3,
    "Ti": 4, "Zr": 4, "Hf": 4, "Sn": 4, "Ge": 4, "Ce": 4, "Ru": 4, "Ir": 4,
    "Nb": 5, "Ta": 5, "V": 5, "P": 5, "Sb": 5,
    "Mo": 6, "W": 6
}

# Extensive library of laboratory-grade precursors for electroceramics
PRECURSOR_DATABASE = {
    "Fe": {
        "Iron(III) Nitrate Nonahydrate (Fe(NO3)3 · 9H2O)": {"mw": 404.00, "phase": "solid", "hygroscopic": True},
        "Iron(III) Chloride (anhydrous, FeCl3)": {"mw": 162.20, "phase": "solid", "hygroscopic": True},
        "Iron(II) Acetate (Fe(CH3COO)2)": {"mw": 173.93, "phase": "solid", "hygroscopic": False}
    },
    "Co": {
        "Cobalt(II) Nitrate Hexahydrate (Co(NO3)2 · 6H2O)": {"mw": 291.03, "phase": "solid", "hygroscopic": True},
        "Cobalt(II) Acetate Tetrahydrate (Co(CH3COO)2 · 4H2O)": {"mw": 249.08, "phase": "solid", "hygroscopic": False}
    },
    "La": {
        "Lanthanum(III) Nitrate Hexahydrate (La(NO3)3 · 6H2O)": {"mw": 433.01, "phase": "solid", "hygroscopic": True},
        "Lanthanum(III) Nitrate (anhydrous, La(NO3)3)": {"mw": 324.91, "phase": "solid", "hygroscopic": True},
        "Lanthanum(III) Acetate Hydrate (La(CH3COO)3 · xH2O)": {"mw": 343.04, "phase": "solid", "hygroscopic": False}
    },
    "Sr": {
        "Strontium Nitrate (Sr(NO3)2, anhydrous)": {"mw": 211.63, "phase": "solid", "hygroscopic": False},
        "Strontium Acetate (Sr(CH3COO)2)": {"mw": 205.71, "phase": "solid", "hygroscopic": False}
    },
    "Ba": {
        "Barium Nitrate (Ba(NO3)2, anhydrous)": {"mw": 261.34, "phase": "solid", "hygroscopic": False},
        "Barium Acetate (Ba(CH3COO)2)": {"mw": 255.42, "phase": "solid", "hygroscopic": False}
    },
    "Ni": {
        "Nickel(II) Nitrate Hexahydrate (Ni(NO3)2 · 6H2O)": {"mw": 290.79, "phase": "solid", "hygroscopic": True},
        "Nickel(II) Acetate Tetrahydrate (Ni(CH3COO)2 · 4H2O)": {"mw": 248.86, "phase": "solid", "hygroscopic": False}
    },
    "Mn": {
        "Manganese(II) Nitrate Tetrahydrate (Mn(NO3)2 · 4H2O)": {"mw": 251.01, "phase": "solid", "hygroscopic": True},
        "Manganese(II) Acetate Tetrahydrate (Mn(CH3COO)2 · 4H2O)": {"mw": 245.09, "phase": "solid", "hygroscopic": False}
    },
    "Ti": {
        "Titanium(IV) Butoxide (Ti(OBu)4, 97%)": {"mw": 340.35, "phase": "liquid", "density": 1.000, "hydrolysis_sensitive": True},
        "Titanium(IV) Isopropoxide (Ti(OiPr)4, 97%)": {"mw": 284.22, "phase": "liquid", "density": 0.960, "hydrolysis_sensitive": True}
    },
    "Zr": {
        "Zirconium(IV) Propoxide (70 wt% in 1-propanol)": {"mw": 327.57 / 0.70, "phase": "liquid", "density": 1.044, "hydrolysis_sensitive": True},
        "Zirconyl Nitrate Hydrate (ZrO(NO3)2 · xH2O)": {"mw": 231.23, "phase": "solid", "hygroscopic": True}
    },
    "Ce": {
        "Cerium(III) Nitrate Hexahydrate (Ce(NO3)3 · 6H2O)": {"mw": 434.22, "phase": "solid", "hygroscopic": True},
        "Ammonium Cerium(IV) Nitrate ((NH4)2Ce(NO3)6)": {"mw": 548.22, "phase": "solid", "hygroscopic": False}
    },
    "Gd": {
        "Gadolinium(III) Nitrate Hexahydrate (Gd(NO3)3 · 6H2O)": {"mw": 451.36, "phase": "solid", "hygroscopic": True}
    },
    "Sm": {
        "Samarium(III) Nitrate Hexahydrate (Sm(NO3)3 · 6H2O)": {"mw": 444.47, "phase": "solid", "hygroscopic": True}
    },
    "Pr": {
        "Praseodymium(III) Nitrate Hexahydrate (Pr(NO3)3 · 6H2O)": {"mw": 435.01, "phase": "solid", "hygroscopic": True}
    },
    "Y": {
        "Yttrium(III) Nitrate Hexahydrate (Y(NO3)3 · 6H2O)": {"mw": 383.01, "phase": "solid", "hygroscopic": True}
    },
    "Li": {
        "Lithium Nitrate (LiNO3, anhydrous)": {"mw": 68.946, "phase": "solid", "hygroscopic": True},
        "Lithium Acetate Dihydrate (CH3COOLi · 2H2O)": {"mw": 102.02, "phase": "solid", "hygroscopic": False},
        "Lithium Carbonate (Li2CO3)": {"mw": 73.891, "phase": "solid", "hygroscopic": False, "stoich_factor": 0.5}
    },
    "Al": {
        "Aluminum Nitrate Nonahydrate (Al(NO3)3 · 9H2O)": {"mw": 375.13, "phase": "solid", "hygroscopic": True},
        "Aluminum Isopropoxide (Al(OiPr)3)": {"mw": 204.24, "phase": "solid", "hydrolysis_sensitive": True}
    },
    "P": {
        "Ammonium Dihydrogen Phosphate ((NH4)H2PO4)": {"mw": 115.03, "phase": "solid", "hygroscopic": False}
    },
    "Cu": {
        "Copper(II) Nitrate Trihydrate (Cu(NO3)2 · 3H2O)": {"mw": 241.60, "phase": "solid", "hygroscopic": True}
    },
    "Zn": {
        "Zinc Nitrate Hexahydrate (Zn(NO3)2 · 6H2O)": {"mw": 297.49, "phase": "solid", "hygroscopic": True}
    }
}

ORGANIC_REAGENTS = {
    "Citric Acid (anhydrous)": {"mw": 192.124},
    "Citric Acid (monohydrate)": {"mw": 210.140},
    "Ethylene Glycol": {"mw": 62.070, "density": 1.1132}
}

# ---------------------------------------------------------
# CHEMICAL PARSING & LOGIC
# ---------------------------------------------------------
def parse_chemical_formula(formula_str):
    cleaned = re.sub(r'[\s_\{\}\$]', '', formula_str)
    pattern = r'([A-Z][a-z]?)([0-9]*\.?[0-9]*)'
    matches = re.findall(pattern, cleaned)
    
    parsed = {}
    for elem, coeff in matches:
        if not elem:
            continue
        val = float(coeff) if coeff else 1.0
        parsed[elem] = parsed.get(elem, 0.0) + val
    return parsed

def compute_formula_weight(composition_dict):
    total_mw = 0.0
    for elem, coeff in composition_dict.items():
        if elem not in ATOMIC_WEIGHTS:
            raise ValueError(f"Unknown element symbol '{elem}'.")
        total_mw += coeff * ATOMIC_WEIGHTS[elem]
    return total_mw

def get_recommended_ratios(composition_dict):
    total_cations = 0.0
    weighted_valency = 0.0
    
    for elem, coeff in composition_dict.items():
        if elem != "O":
            total_cations += coeff
            val = STANDARD_VALENCIES.get(elem, 3)
            weighted_valency += coeff * val
            
    avg_charge = weighted_valency / total_cations if total_cations > 0 else 3.0
    rec_ca_ratio = max(2.5, 1.0 + 0.5 * avg_charge)
    rec_eg_ratio = 6.0 if avg_charge <= 2.5 else 8.0
    return avg_charge, rec_ca_ratio, rec_eg_ratio

# ---------------------------------------------------------
# INTERACTIVE APPLICATION UI
# ---------------------------------------------------------
col1, col2 = st.columns([1.2, 1.8], gap="large")

with col1:
    st.subheader("1. Material System & Presets")
    preset_choice = st.selectbox(
        "Choose Common Electrochemical Benchmark or Custom:",
        [
            "LaCo0.8Fe0.2O3 (SOFC / Electrocatalyst Perovskite)",
            "La0.6Sr0.4Co0.2Fe0.8O3 (LSCF Cathode)",
            "La0.8Sr0.2MnO3 (LSM SOFC Cathode)",
            "Ce0.9Gd0.1O1.95 (GDC Electrolyte)",
            "BaZr0.8Y0.2O3 (BZY Proton Conductor)",
            "Li0.40La0.60TiO3 (LLTO Perovskite Battery Electrolyte)",
            "Li7La3Zr2O12 (LLZO Garnet Battery Electrolyte)",
            "LiNi0.5Mn1.5O4 (High-Voltage Battery Spinel)",
            "Custom Composition"
        ]
    )
    
    preset_map = {
        "LaCo0.8Fe0.2O3": ("LaCo0.8Fe0.2O3", 950),
        "LSCF": ("La0.6Sr0.4Co0.2Fe0.8O3", 1000),
        "LSM": ("La0.8Sr0.2MnO3", 1000),
        "GDC": ("Ce0.9Gd0.1O1.95", 800),
        "BZY": ("BaZr0.8Y0.2O3", 1100),
        "LLTO": ("Li0.40La0.60TiO3", 850),
        "LLZO": ("Li7La3Zr2O12", 1100),
        "LNMO": ("LiNi0.5Mn1.5O4", 800),
    }

    selected_formula = "LaCo0.8Fe0.2O3"
    selected_temp = 950
    for key, (form, temp) in preset_map.items():
        if key in preset_choice:
            selected_formula = form
            selected_temp = temp
            break

    formula_str = st.text_input("Chemical Formula:", value=selected_formula)
    batch_mass_g = st.number_input("Target Sintered Powder Mass (g):", min_value=0.1, max_value=1000.0, value=10.0, step=0.5)
    sinter_temp_c = st.slider("Final Sintering / Calcination Temp (°C):", 500, 1400, selected_temp, 25)

    # Volatile cation loss rule
    has_li = "Li" in parse_chemical_formula(formula_str)
    if has_li:
        if sinter_temp_c <= 750:
            auto_li_excess = 0.0
        elif sinter_temp_c <= 900:
            auto_li_excess = 5.0 * ((sinter_temp_c - 750) / 150.0)
        else:
            auto_li_excess = 5.0 + 10.0 * ((sinter_temp_c - 900) / 300.0)
        li_excess_pct = st.number_input("Volatile Li Sublimation Compensation (%):", 0.0, 50.0, auto_li_excess, 1.0)
    else:
        li_excess_pct = 0.0

with col2:
    st.subheader("2. Coordination & Precursor Matrix")
    
    try:
        comp = parse_chemical_formula(formula_str)
        mw_oxide = compute_formula_weight(comp)
        avg_charge, rec_ca, rec_eg = get_recommended_ratios(comp)
        
        m1, m2, m3 = st.columns(3)
        m1.metric("Unit MW", f"{mw_oxide:.3f} g/mol")
        m2.metric("Mean Cation Charge", f"+{avg_charge:.2f}")
        m3.metric("Nominal Moles", f"{(batch_mass_g / mw_oxide):.5f} mol")
        
    except Exception as err:
        st.error(f"Formula Error: {err}")
        st.stop()

    st.markdown("#### Chelation & Polymerization Ratios")
    c_che1, c_che2 = st.columns(2)
    with c_che1:
        ca_type = st.selectbox("Citric Acid Form:", ["Citric Acid (anhydrous)", "Citric Acid (monohydrate)"])
        ca_ratio = st.slider("Citric Acid : Total Cation Molar Ratio:", 1.5, 6.0, float(np.round(rec_ca, 1)), 0.1)
    with c_che2:
        eg_ratio = st.slider("Ethylene Glycol : Citric Acid Molar Ratio:", 2.0, 14.0, float(rec_eg), 0.5)

    st.markdown("#### Precursor Sourcing (With Dynamic Fallback for Any Element)")
    metal_cations = [e for e in comp.keys() if e != "O"]
    
    selected_precursors = {}
    
    cols = st.columns(min(len(metal_cations), 4))
    for idx, elem in enumerate(metal_cations):
        col_target = cols[idx % len(cols)]
        with col_target:
            if elem in PRECURSOR_DATABASE:
                reagent_options = list(PRECURSOR_DATABASE[elem].keys())
                chosen = st.selectbox(f"{elem} Reagent:", reagent_options, key=f"sel_{elem}")
                selected_precursors[elem] = PRECURSOR_DATABASE[elem][chosen]
                selected_precursors[elem]["name"] = chosen
            else:
                st.warning(f"Using Universal Fallback for {elem}")
                default_mw = ATOMIC_WEIGHTS.get(elem, 100.0) + (3 * 62.00) # Default to nitrate
                custom_name = st.text_input(f"{elem} Reagent Name:", value=f"{elem}(NO3)3 hydrate", key=f"name_{elem}")
                custom_mw = st.number_input(f"{elem} MW (g/mol):", value=default_mw, key=f"mw_{elem}")
                is_liquid = st.checkbox(f"{elem} is Liquid?", value=False, key=f"liq_{elem}")
                density = st.number_input(f"{elem} Density (g/mL):", value=1.00, key=f"den_{elem}") if is_liquid else 1.00
                selected_precursors[elem] = {
                    "name": custom_name,
                    "mw": custom_mw,
                    "phase": "liquid" if is_liquid else "solid",
                    "density": density,
                    "hygroscopic": True
                }

# ---------------------------------------------------------
# GENERATE DISPENSING SHEET & SOP
# ---------------------------------------------------------
st.markdown("---")
if st.button("Generate Experimental Bench SOP", type="primary"):
    n_target = batch_mass_g / mw_oxide
    
    elem_moles = {}
    for elem, coeff in comp.items():
        if elem == "O":
            continue
        if elem == "Li":
            elem_moles[elem] = n_target * coeff * (1.0 + (li_excess_pct / 100.0))
        else:
            elem_moles[elem] = n_target * coeff

    total_cation_moles = sum(elem_moles.values())
    
    precursor_table = []
    moisture_warnings = []
    
    for elem, moles in elem_moles.items():
        reagent = selected_precursors[elem]
        stoich_fac = reagent.get("stoich_factor", 1.0)
        reagent_moles = moles / stoich_fac
        mass_g = reagent_moles * reagent["mw"]
        
        if reagent["phase"] == "liquid":
            vol_ml = mass_g / reagent["density"]
            dispense_action = f"**{vol_ml:.2f} mL** micropipette ({mass_g:.3f} g)"
        else:
            dispense_action = f"**{mass_g:.3f} g** balance"
            
        if reagent.get("hydrolysis_sensitive"):
            moisture_warnings.append(f"**{elem} ({reagent['name']})**: Rapid hydrolysis! Do not expose to water. Dissolve into acidic glycol sol.")
            
        precursor_table.append({
            "Element": elem,
            "Reagent": reagent["name"],
            "Moles": f"{reagent_moles:.5f} mol",
            "Dispense Quantity": dispense_action,
            "Handling Note": "Hygroscopic ⚠️" if reagent.get("hygroscopic") else "Standard"
        })

    moles_ca = total_cation_moles * ca_ratio
    mass_ca = moles_ca * ORGANIC_REAGENTS[ca_type]["mw"]
    
    moles_eg = moles_ca * eg_ratio
    mass_eg = moles_eg * ORGANIC_REAGENTS["Ethylene Glycol"]["mw"]
    vol_eg = mass_eg / ORGANIC_REAGENTS["Ethylene Glycol"]["density"]
    est_nh3_vol = (moles_ca * 1.2 * 17.031) / (0.90 * 0.28)

    organics_table = [
        {"Reagent": ca_type, "Moles": f"{moles_ca:.4f} mol", "Dispense Quantity": f"**{mass_ca:.2f} g**", "Role": "Polydentate Chelation"},
        {"Reagent": "Ethylene Glycol", "Moles": f"{moles_eg:.4f} mol", "Dispense Quantity": f"**{vol_eg:.1f} mL** ({mass_eg:.2f} g)", "Role": "Polyester Polymerization"},
        {"Reagent": "Aqueous Ammonia (NH4OH, ~28%)", "Moles": f"~{moles_ca*1.2:.3f} mol", "Dispense Quantity": f"**~{est_nh3_vol:.1f} mL** (to pH 7.0)", "Role": "Deprotonation / Buffering"}
    ]

    st.subheader("📋 Experimental Dispensing Sheet")
    col_t1, col_t2 = st.columns([1.3, 1.0])
    with col_t1:
        st.markdown("#### Inorganic Metal Precursors")
        st.dataframe(pd.DataFrame(precursor_table), use_container_width=True, hide_index=True)
    with col_t2:
        st.markdown("#### Organic Chelation & Buffer Matrix")
        st.dataframe(pd.DataFrame(organics_table), use_container_width=True, hide_index=True)

    if moisture_warnings:
        st.warning("⚠️ **Moisture Sensitivity Alerts:**\n\n" + "\n\n".join(moisture_warnings))

    st.subheader("🧪 Standard Operating Procedure (Step-by-Step Bench Protocol)")
    st.markdown(f"""
    1. **Ligand Solvation:**
       - Measure **{vol_eg:.1f} mL** of Ethylene Glycol. Divide roughly 50:50 into Beaker A and Beaker B.
       - Dissolve **{mass_ca:.2f} g** of {ca_type} in Beaker A at **50–60°C** under stirring until clear.
    2. **Metal Salt Chelation:**
       - Dissolve the nitrate / acetate salts in Beaker B.
       - Combine Beaker A and B. Stir at 50°C for 30 minutes to form stable metal-citrate complexes. *(If using moisture-sensitive alkoxides, inject them slowly into the clear solution only after it has cooled below 50°C)*.
    3. **Deprotonation & pH Adjustment:**
       - Add aqueous $\\text{{NH}}_4\\text{{OH}}$ dropwise until reaching **pH = 7.0 ± 0.2**. This deprotonates the carboxylic acid groups, locking cations into the spatial chelation complex.
    4. **Polyesterification & Char Formation:**
       - Heat at **140–150°C** to drive off water and promote ester condensation, yielding a viscous, amber polyester gel.
       - Pyrolyze in a crucible at **200–250°C for 30–60 minutes** until it expands into a black solid precursor foam.
    5. **Calcination & Crystallization:**
       - Heat at **450–500°C for 5 hours** to burn out all residual carbonaceous matter.
       - Grind the precursor in an agate mortar or planetary ball mill.
       - Sinter at **{sinter_temp_c}°C for 8–24 hours** in air (or controlled atmosphere) to obtain single-phase crystalline **{formula_str}**.
    """)
