import streamlit as st
import pandas as pd
import numpy as np
import re

# ---------------------------------------------------------
# PAGE SETUP & SCIENTIFIC STYLING
# ---------------------------------------------------------
st.set_page_config(
    page_title="PechiniLab | Electrochemistry Sol-Gel Suite",
    page_icon="⚗️",
    layout="wide"
)

st.title("⚗️ PechiniLab: Precision Sol-Gel Engine for Electroceramics")
st.markdown(
    "A rigorous, thermodynamic- and stoichiometry-driven planning suite for **polymeric precursor (Pechini) "
    "synthesis** of solid-state battery electrolytes and cathode oxides."
)

# ---------------------------------------------------------
# PECHINI & SOLID-STATE KNOWLEDGE BASE
# ---------------------------------------------------------
ATOMIC_WEIGHTS = {
    "H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998,
    "Li": 6.941, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085,
    "P": 30.974, "S": 32.065, "Cl": 35.453, "K": 39.098, "Ca": 40.078,
    "Sc": 44.956, "Ti": 47.867, "V": 50.942, "Cr": 51.996, "Mn": 54.938,
    "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546, "Zn": 65.380,
    "Ga": 69.723, "Ge": 72.630, "Sr": 87.620, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.950, "Sn": 118.710, "Ba": 137.327, "La": 138.905,
    "Ce": 140.116, "Pr": 140.908, "Nd": 144.242, "Sm": 150.360, "Eu": 151.964,
    "Gd": 157.250, "Ta": 180.948, "W": 183.840, "Bi": 208.980
}

CATION_VALENCIES = {
    "Li": 1, "Na": 1, "K": 1, "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "Zn": 2,
    "Al": 3, "Sc": 3, "Y": 3, "La": 3, "Ga": 3, "Fe": 3, "Cr": 3, "Ce": 3,
    "Ti": 4, "Zr": 4, "Sn": 4, "Ge": 4, "Mn": 4, "Nb": 5, "Ta": 5, "P": 5,
    "Ni": 2, "Co": 3
}

# Precursor Library: molecular weight, physical phase, density (g/mL), moisture sensitivity
PRECURSOR_DATABASE = {
    "Li": {
        "Lithium Nitrate (LiNO3, anhydrous)": {"formula": "LiNO3", "mw": 68.946, "phase": "solid", "hygroscopic": True},
        "Lithium Acetate Dihydrate (CH3COOLi · 2H2O)": {"formula": "C2H7LiO4", "mw": 102.02, "phase": "solid", "hygroscopic": False},
        "Lithium Carbonate (Li2CO3)": {"formula": "Li2CO3", "mw": 73.891, "phase": "solid", "hygroscopic": False, "stoich_factor": 0.5}
    },
    "La": {
        "Lanthanum(III) Nitrate Hexahydrate (La(NO3)3 · 6H2O)": {"formula": "La(NO3)3·6H2O", "mw": 433.01, "phase": "solid", "hygroscopic": True},
        "Lanthanum(III) Nitrate (anhydrous)": {"formula": "La(NO3)3", "mw": 324.91, "phase": "solid", "hygroscopic": True},
        "Lanthanum(III) Acetate Hydrate (La(CH3COO)3 · xH2O)": {"formula": "La(C2H3O2)3·1.5H2O", "mw": 343.04, "phase": "solid", "hygroscopic": False}
    },
    "Ti": {
        "Titanium(IV) Butoxide (Ti(OBu)4, 97%)": {"formula": "Ti(OCH2CH2CH2CH3)4", "mw": 340.35, "phase": "liquid", "density": 1.000, "hydrolysis_sensitive": True},
        "Titanium(IV) Isopropoxide (Ti(OiPr)4, 97%)": {"formula": "Ti(OCH(CH3)2)4", "mw": 284.22, "phase": "liquid", "density": 0.960, "hydrolysis_sensitive": True}
    },
    "Zr": {
        "Zirconium(IV) Propoxide (70 wt% in 1-propanol)": {"formula": "Zr(OC3H7)4", "mw": 327.57 / 0.70, "phase": "liquid", "density": 1.044, "hydrolysis_sensitive": True},
        "Zirconyl Nitrate Hydrate (ZrO(NO3)2 · xH2O)": {"formula": "ZrO(NO3)2·xH2O", "mw": 231.23, "phase": "solid", "hygroscopic": True}
    },
    "Al": {
        "Aluminum Nitrate Nonahydrate (Al(NO3)3 · 9H2O)": {"formula": "Al(NO3)3·9H2O", "mw": 375.13, "phase": "solid", "hygroscopic": True},
        "Aluminum Isopropoxide (Al(OiPr)3)": {"formula": "Al(C3H7O)3", "mw": 204.24, "phase": "solid", "hydrolysis_sensitive": True}
    },
    "P": {
        "Ammonium Dihydrogen Phosphate ((NH4)H2PO4)": {"formula": "(NH4)H2PO4", "mw": 115.03, "phase": "solid", "hygroscopic": False}
    },
    "Sr": {
        "Strontium Nitrate (Sr(NO3)2, anhydrous)": {"formula": "Sr(NO3)2", "mw": 211.63, "phase": "solid", "hygroscopic": False}
    },
    "Ba": {
        "Barium Nitrate (Ba(NO3)2, anhydrous)": {"formula": "Ba(NO3)2", "mw": 261.34, "phase": "solid", "hygroscopic": False}
    },
    "Ni": {
        "Nickel(II) Nitrate Hexahydrate (Ni(NO3)2 · 6H2O)": {"formula": "Ni(NO3)2·6H2O", "mw": 290.79, "phase": "solid", "hygroscopic": True}
    },
    "Mn": {
        "Manganese(II) Nitrate Tetrahydrate (Mn(NO3)2 · 4H2O)": {"formula": "Mn(NO3)2·4H2O", "mw": 251.01, "phase": "solid", "hygroscopic": True}
    },
    "Co": {
        "Cobalt(II) Nitrate Hexahydrate (Co(NO3)2 · 6H2O)": {"formula": "Co(NO3)2·6H2O", "mw": 291.03, "phase": "solid", "hygroscopic": True}
    }
}

ORGANIC_REAGENTS = {
    "Citric Acid (anhydrous)": {"mw": 192.124, "equiv_COOH": 3},
    "Citric Acid (monohydrate)": {"mw": 210.140, "equiv_COOH": 3},
    "Ethylene Glycol": {"mw": 62.070, "density": 1.1132, "equiv_OH": 2}
}

# ---------------------------------------------------------
# CHEMICAL LOGIC ENGINE
# ---------------------------------------------------------
def parse_chemical_formula(formula_str):
    """
    Parses chemical formulas with decimal coefficients (e.g., Li0.4La0.6TiO3, Li7La3Zr2O12).
    Returns dictionary of {Element: coefficient}.
    """
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
    """Calculates molecular weight from elemental dictionary."""
    total_mw = 0.0
    for elem, coeff in composition_dict.items():
        if elem not in ATOMIC_WEIGHTS:
            raise ValueError(f"Unknown element '{elem}' in formula.")
        total_mw += coeff * ATOMIC_WEIGHTS[elem]
    return total_mw

def get_recommended_ratios(composition_dict):
    """
    Calculates thermodynamically stable chelation (CA:Metal) and
    polyesterification (EG:CA) ratios based on average cation charge.
    """
    total_cations = 0.0
    weighted_valency = 0.0
    
    for elem, coeff in composition_dict.items():
        if elem != "O":
            total_cations += coeff
            val = CATION_VALENCIES.get(elem, 2)
            weighted_valency += coeff * val
            
    avg_charge = weighted_valency / total_cations if total_cations > 0 else 2.0
    # Higher average cation charge demands higher ligand density to prevent hydrolysis
    rec_ca_ratio = max(2.5, 1.0 + 0.5 * avg_charge)
    rec_eg_ratio = 6.0 if avg_charge <= 2.5 else 8.0
    return avg_charge, rec_ca_ratio, rec_eg_ratio

# ---------------------------------------------------------
# USER INTERFACE: SYNTHESIS PLANNER
# ---------------------------------------------------------
col_left, col_right = st.columns([1.2, 1.8], gap="large")

with col_left:
    st.subheader("1. Composition & Batch Parameters")
    
    preset_selection = st.selectbox(
        "Choose Composition Preset or Custom:",
        [
            "Li0.40La0.60TiO3 (LLTO Perovskite)",
            "Li7La3Zr2O12 (LLZO Garnet)",
            "Li1.3Al0.3Ti1.7(PO4)3 (LATP NASICON)",
            "LiNi0.5Mn1.5O4 (High-Voltage Spinel LNMO)",
            "Custom Target Formulation"
        ]
    )
    
    if "LLTO" in preset_selection:
        default_formula = "Li0.40La0.60TiO3"
        default_temp = 850
    elif "LLZO" in preset_selection:
        default_formula = "Li7La3Zr2O12"
        default_temp = 1100
    elif "LATP" in preset_selection:
        default_formula = "Li1.3Al0.3Ti1.7P3O12"
        default_temp = 900
    elif "LNMO" in preset_selection:
        default_formula = "LiNi0.5Mn1.5O4"
        default_temp = 800
    else:
        default_formula = "Li0.33La0.55TiO3"
        default_temp = 850
        
    formula_input = st.text_input("Chemical Formula:", value=default_formula)
    batch_mass_g = st.number_input("Target Sintered Ceramic Yield (g):", min_value=0.1, max_value=1000.0, value=10.0, step=1.0)
    
    st.markdown("#### High-Temperature Kinetics")
    sinter_temp_c = st.slider("Final Calcination / Sintering Temp (°C):", min_value=600, max_value=1300, value=default_temp, step=25)
    
    # Calculate automated lithium loss compensation
    if sinter_temp_c <= 750:
        auto_li_excess = 0.0
    elif sinter_temp_c <= 900:
        auto_li_excess = 5.0 * ((sinter_temp_c - 750) / 150.0)
    else:
        auto_li_excess = 5.0 + 10.0 * ((sinter_temp_c - 900) / 300.0)
        
    override_excess = st.checkbox("Manually override volatile Li compensation", value=False)
    if override_excess:
        li_excess_pct = st.number_input("Lithium Volatilization Excess (%):", min_value=0.0, max_value=50.0, value=10.0, step=1.0)
    else:
        li_excess_pct = auto_li_excess
        st.info(f"Knudsen-Langmuir Li2O sublimation compensation: **+{li_excess_pct:.1f}% excess Li** at {sinter_temp_c}°C.")

with col_right:
    st.subheader("2. Complexation & Polymerization Dynamics")
    
    try:
        composition = parse_chemical_formula(formula_input)
        mw_target = compute_formula_weight(composition)
        avg_charge, rec_ca, rec_eg = get_recommended_ratios(composition)
        
        c_stat1, c_stat2, c_stat3 = st.columns(3)
        c_stat1.metric("Formula Weight", f"{mw_target:.3f} g/mol")
        c_stat2.metric("Mean Cation Charge (z̄)", f"+{avg_charge:.2f}")
        c_stat3.metric("Nominal Moles Target", f"{(batch_mass_g / mw_target):.5f} mol")
        
    except Exception as e:
        st.error(f"Error parsing chemical formula: {e}")
        st.stop()

    st.markdown("#### Chelation & Transesterification Network")
    col_che1, col_che2 = st.columns(2)
    with col_che1:
        ca_type = st.selectbox("Citric Acid Hydration:", ["Citric Acid (anhydrous)", "Citric Acid (monohydrate)"])
        ca_ratio = st.slider("Citric Acid : Total Cation Molar Ratio:", min_value=1.5, max_value=5.0, value=float(np.round(rec_ca, 1)), step=0.1,
                             help="Controls chelation stability. Ratios < 2.0 risk precipitation during pH buffering.")
    with col_che2:
        eg_ratio = st.slider("Ethylene Glycol : Citric Acid Molar Ratio:", min_value=2.0, max_value=12.0, value=float(rec_eg), step=0.5,
                             help="Controls polyester crosslink density. Ratios < 3 fail to crosslink; ratios > 10 create uncontrolled combustion exotherms.")

    st.markdown("#### Precursor Reagent Selection")
    selected_reagents = {}
    missing_elements = []
    
    metal_cations = [elem for elem in composition.keys() if elem != "O"]
    precursor_cols = st.columns(len(metal_cations))
    
    for i, elem in enumerate(metal_cations):
        with precursor_cols[i]:
            if elem in PRECURSOR_DATABASE:
                reagent_options = list(PRECURSOR_DATABASE[elem].keys())
                selected_reagents[elem] = st.selectbox(f"{elem} Source:", reagent_options)
            else:
                missing_elements.append(elem)
                st.error(f"No precursor profile for: {elem}")
                
    if missing_elements:
        st.stop()

# ---------------------------------------------------------
# CALCULATION & RECIPE GENERATION
# ---------------------------------------------------------
st.markdown("---")
if st.button("Generate Experimental Bench SOP", type="primary"):
    n_target = batch_mass_g / mw_target
    
    # 1. Calculate elemental moles
    elem_moles = {}
    for elem, coeff in composition.items():
        if elem == "O":
            continue
        if elem == "Li":
            elem_moles[elem] = n_target * coeff * (1.0 + (li_excess_pct / 100.0))
        else:
            elem_moles[elem] = n_target * coeff
            
    total_cation_moles = sum(elem_moles.values())
    
    # 2. Build Precursor Dispensing Table
    precursor_results = []
    hydrolysis_warnings = []
    
    for elem, moles in elem_moles.items():
        chosen_name = selected_reagents[elem]
        spec = PRECURSOR_DATABASE[elem][chosen_name]
        
        # Adjust for stoichiometry factors (e.g. Li2CO3 has 2 Li per mole)
        stoich_fac = spec.get("stoich_factor", 1.0)
        reagent_moles = moles / stoich_fac
        mass_reagent_g = reagent_moles * spec["mw"]
        
        if spec["phase"] == "liquid":
            vol_mL = mass_reagent_g / spec["density"]
            dispense_action = f"**{vol_mL:.2f} mL** via micropipette ({mass_reagent_g:.3f} g)"
        else:
            dispense_action = f"**{mass_reagent_g:.3f} g** on analytical balance"
            
        if spec.get("hydrolysis_sensitive", False):
            hydrolysis_warnings.append(f"**{elem} ({chosen_name})**: Reacts violently with moisture. Inject slowly into EG/CA sol only after cooling below 60°C.")
            
        precursor_results.append({
            "Cation": elem,
            "Reagent Selected": chosen_name,
            "Molar Requirement": f"{reagent_moles:.5f} mol",
            "Dispense Quantity": dispense_action,
            "Hygroscopic / Volatile Flag": "Hygroscopic ⚠️" if spec.get("hygroscopic") else ("Volatile Sublimation ⚠️" if elem == "Li" and li_excess_pct > 0 else "Stable")
        })

    # 3. Organics & Buffer Calculations
    moles_ca = total_cation_moles * ca_ratio
    mass_ca = moles_ca * ORGANIC_REAGENTS[ca_type]["mw"]
    
    moles_eg = moles_ca * eg_ratio
    mass_eg = moles_eg * ORGANIC_REAGENTS["Ethylene Glycol"]["mw"]
    vol_eg = mass_eg / ORGANIC_REAGENTS["Ethylene Glycol"]["density"]
    
    # Ammonia buffer estimate (typically 1.0 to 1.5 moles NH3 per mole CA to reach pH 7)
    est_nh3_vol_ml = (moles_ca * 1.2 * 17.031) / (0.90 * 0.28) # 28% aqueous NH3
    
    organics_results = [
        {"Component": ca_type, "Moles": f"{moles_ca:.4f} mol", "Dispense Quantity": f"**{mass_ca:.2f} g**", "Function": "Polydentate Chelation Matrix"},
        {"Component": "Ethylene Glycol", "Moles": f"{moles_eg:.4f} mol", "Dispense Quantity": f"**{vol_eg:.1f} mL** ({mass_eg:.2f} g)", "Function": "Transesterification Crosslinker"},
        {"Component": "Aqueous Ammonia (NH4OH, ~28%)", "Moles": f"~{moles_ca*1.2:.3f} mol", "Dispense Quantity": f"**~{est_nh3_vol_ml:.1f} mL** (dropwise to pH 7.0)", "Function": "Carboxylic Deprotonation Buffer"}
    ]

    # ----------------- DISPLAY TABLES -----------------
    st.subheader("📋 Experimental Dispensing Sheet")
    
    tab_a, tab_b = st.columns([1.3, 1.0], gap="medium")
    with tab_a:
        st.markdown("#### Inorganic Precursor Salts & Alkoxides")
        st.dataframe(pd.DataFrame(precursor_results), use_container_width=True, hide_index=True)
    with tab_b:
        st.markdown("#### Organic Complexation & Buffer Matrix")
        st.dataframe(pd.DataFrame(organics_results), use_container_width=True, hide_index=True)

    if hydrolysis_warnings:
        st.warning("⚠️ **Moisture Sensitivity Alerts:**\n\n" + "\n\n".join(hydrolysis_warnings))

    # ----------------- STEP-BY-STEP SOP -----------------
    st.subheader("🧪 Standard Operating Procedure (Step-by-Step Bench Protocol)")
    
    st.markdown(f"""
    1. **Pre-Dissolution & Ligand Solvation:**
       - Divide the calculated **{vol_eg:.1f} mL** of Ethylene Glycol into two separate beakers (Beaker A: ~50%, Beaker B: ~50%).
       - Add **{mass_ca:.2f} g** of {ca_type} into Beaker A. Heat to **50–60°C** under continuous magnetic stirring until completely transparent.
    
    2. **Cation Chelation Sequence:**
       - Dissolve the non-hydrolyzing metal salts into Beaker B[cite: 1].
       - Combine Beaker A and Beaker B under vigorous stirring. Allow the temperature to equilibrate to ~40–50°C[cite: 1].
       - *Alkoxide Injection (if applicable):* Using a dry micropipette, inject liquid alkoxides slowly into the vortex. The solution must remain clear with zero opalescence[cite: 1].
    
    3. **Deprotonation & pH Neutralization:**
       - Insert a calibrated pH electrode into the room-temperature sol[cite: 1].
       - Add aqueous $\\text{{NH}}_4\\text{{OH}}$ dropwise with continuous stirring until reaching **pH = 7.0 ± 0.2**[cite: 1]. *(Note: Neutralization deprotonates $-\\text{{COOH}}$ groups into carboxylates, completing metal-ligand chelation and preventing cation phase separation)[cite: 1].*
    
    4. **Polyesterification & Char Formation:**
       - Heat the neutralized solution to **140–150°C** on a hotplate. Excess water and free glycol will evaporate, driving the transesterification reaction to yield a transparent, viscous amber resin[cite: 1].
       - Transfer the gel into an alumina crucible and heat in a box furnace to **200–230°C for 30 minutes** to expand into a porous solid black foam[cite: 1].
    
    5. **Pyrolysis & Sintering Thermal Profile:**
       - **Organic Burnout:** Ramp at $5^\\circ\\text{{C/min}}$ to **$450^\\circ\\text{{C}}$ for 5 hours** to achieve complete decarbonation and eliminate residual organic matrices[cite: 1].
       - **Intermediate Milling:** Grind the resulting white/grey powder in an agate mortar or planetary ball mill to break soft agglomerates[cite: 1].
       - **Final Crystallization:** Sinter at **{sinter_temp_c}°C for 24–28 hours** in air to achieve full phase purity of **{formula_input}**[cite: 1].
    """)
