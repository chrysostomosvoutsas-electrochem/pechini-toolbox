import streamlit as st
import pandas as pd
import numpy as np
import re
from datetime import datetime

# -----------------------------------------------------------------------------
# APPLICATION CONFIGURATION & BENCH-READY THEME
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="PechiniLab | Departmental Sol-Gel Platform",
    page_icon="⚗️",
    layout="wide"
)

st.title("⚗️ PechiniLab: Production Sol-Gel Engine for Electroceramics")
st.caption("Validated for Solid Oxide Fuel Cells/Electrolyzers (SOFC/SOEC), Solid-State Batteries (SSB), and Electrocatalysts (AEM/PEM/OER).")

# -----------------------------------------------------------------------------
# THERMODYNAMIC & PERIODIC TABLE DATABASE
# -----------------------------------------------------------------------------
ATOMIC_WEIGHTS = {
    "H": 1.008, "Li": 6.941, "Be": 9.012, "B": 10.811, "C": 12.011, "N": 14.007,
    "O": 15.999, "F": 18.998, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085,
    "P": 30.974, "S": 32.065, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867,
    "V": 50.942, "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693,
    "Cu": 63.546, "Zn": 65.380, "Ga": 69.723, "Ge": 72.630, "Sr": 87.620, "Y": 88.906,
    "Zr": 91.224, "Nb": 92.906, "Mo": 95.950, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42,
    "Ag": 107.87, "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Ba": 137.33, "La": 138.91,
    "Ce": 140.12, "Pr": 140.91, "Nd": 144.24, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25,
    "Tb": 158.93, "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05,
    "Lu": 174.97, "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Ir": 192.22, "Pt": 195.08,
    "Pb": 207.20, "Bi": 208.98
}

CATION_SPECS = {
    "Li": {"valency": 1, "log_k1": 0.6, "alkoxide": False, "volatile": True},
    "Na": {"valency": 1, "log_k1": 0.5, "alkoxide": False, "volatile": True},
    "K":  {"valency": 1, "log_k1": 0.4, "alkoxide": False, "volatile": True},
    "Mg": {"valency": 2, "log_k1": 3.4, "alkoxide": False, "volatile": False},
    "Ca": {"valency": 2, "log_k1": 4.8, "alkoxide": False, "volatile": False},
    "Sr": {"valency": 2, "log_k1": 3.2, "alkoxide": False, "volatile": False},
    "Ba": {"valency": 2, "log_k1": 2.8, "alkoxide": False, "volatile": False, "carbonate_trap": True},
    "Co": {"valency": 2, "log_k1": 4.8, "alkoxide": False, "volatile": False},
    "Ni": {"valency": 2, "log_k1": 5.4, "alkoxide": False, "volatile": False},
    "Cu": {"valency": 2, "log_k1": 5.9, "alkoxide": False, "volatile": False},
    "Zn": {"valency": 2, "log_k1": 4.9, "alkoxide": False, "volatile": False},
    "Fe": {"valency": 3, "log_k1": 11.5, "alkoxide": False, "volatile": False},
    "Mn": {"valency": 3, "log_k1": 3.7, "alkoxide": False, "volatile": False},
    "Al": {"valency": 3, "log_k1": 8.0, "alkoxide": True, "volatile": False},
    "Sc": {"valency": 3, "log_k1": 8.5, "alkoxide": False, "volatile": False},
    "Y":  {"valency": 3, "log_k1": 8.0, "alkoxide": False, "volatile": False},
    "La": {"valency": 3, "log_k1": 8.5, "alkoxide": False, "volatile": False},
    "Ce": {"valency": 3, "log_k1": 8.5, "alkoxide": False, "volatile": False},
    "Pr": {"valency": 3, "log_k1": 8.6, "alkoxide": False, "volatile": False},
    "Nd": {"valency": 3, "log_k1": 8.7, "alkoxide": False, "volatile": False},
    "Sm": {"valency": 3, "log_k1": 8.8, "alkoxide": False, "volatile": False},
    "Gd": {"valency": 3, "log_k1": 8.9, "alkoxide": False, "volatile": False},
    "Ti": {"valency": 4, "log_k1": 18.0, "alkoxide": True, "volatile": False},
    "Zr": {"valency": 4, "log_k1": 16.0, "alkoxide": True, "volatile": False},
    "Sn": {"valency": 4, "log_k1": 14.0, "alkoxide": False, "volatile": False},
    "Nb": {"valency": 5, "log_k1": 15.0, "alkoxide": True, "volatile": False},
    "Ta": {"valency": 5, "log_k1": 15.0, "alkoxide": True, "volatile": False},
    "P":  {"valency": 5, "log_k1": 1.0, "alkoxide": False, "volatile": True}
}

PRECURSOR_LIBRARY = {
    "La": {
        "La(NO3)3 · 6H2O": {"mw": 433.01, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 99.0},
        "La(CH3COO)3 · 1.5H2O": {"mw": 343.04, "phase": "solid", "hygro": False, "type": "acetate", "default_assay": 99.0}
    },
    "Co": {
        "Co(NO3)2 · 6H2O": {"mw": 291.03, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 98.5},
        "Co(CH3COO)2 · 4H2O": {"mw": 249.08, "phase": "solid", "hygro": False, "type": "acetate", "default_assay": 99.0}
    },
    "Fe": {
        "Fe(NO3)3 · 9H2O": {"mw": 404.00, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 98.0},
        "Fe(CH3COO)2": {"mw": 173.93, "phase": "solid", "hygro": False, "type": "acetate", "default_assay": 97.0},
        "FeCl3 (anhydrous)": {"mw": 162.20, "phase": "solid", "hygro": True, "type": "chloride", "default_assay": 97.0}
    },
    "Sr": {
        "Sr(NO3)2 (anhydrous)": {"mw": 211.63, "phase": "solid", "hygro": False, "type": "nitrate", "default_assay": 99.0},
        "Sr(CH3COO)2": {"mw": 205.71, "phase": "solid", "hygro": False, "type": "acetate", "default_assay": 99.0}
    },
    "Ba": {
        "Ba(NO3)2 (anhydrous)": {"mw": 261.34, "phase": "solid", "hygro": False, "type": "nitrate", "default_assay": 99.0},
        "Ba(CH3COO)2": {"mw": 255.42, "phase": "solid", "hygro": False, "type": "acetate", "default_assay": 99.0}
    },
    "Ni": {
        "Ni(NO3)2 · 6H2O": {"mw": 290.79, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 98.5},
        "Ni(CH3COO)2 · 4H2O": {"mw": 248.86, "phase": "solid", "hygro": False, "type": "acetate", "default_assay": 99.0}
    },
    "Mn": {
        "Mn(NO3)2 · 4H2O": {"mw": 251.01, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 98.0},
        "Mn(CH3COO)2 · 4H2O": {"mw": 245.09, "phase": "solid", "hygro": False, "type": "acetate", "default_assay": 99.0}
    },
    "Ti": {
        "Ti(OBu)4 (Titanium(IV) Butoxide, 97%)": {"mw": 340.35, "phase": "liquid", "density": 1.000, "type": "alkoxide", "default_assay": 97.0},
        "Ti(OiPr)4 (Titanium(IV) Isopropoxide, 97%)": {"mw": 284.22, "phase": "liquid", "density": 0.960, "type": "alkoxide", "default_assay": 97.0}
    },
    "Zr": {
        "Zr(OPr)4 (Zirconium Propoxide 70% in 1-propanol)": {"mw": 327.57 / 0.70, "phase": "liquid", "density": 1.044, "type": "alkoxide", "default_assay": 100.0},
        "ZrO(NO3)2 · xH2O (Zirconyl Nitrate)": {"mw": 231.23, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 99.0}
    },
    "Ce": {
        "Ce(NO3)3 · 6H2O": {"mw": 434.22, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 99.0}
    },
    "Gd": {
        "Gd(NO3)3 · 6H2O": {"mw": 451.36, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 99.0}
    },
    "Sm": {
        "Sm(NO3)3 · 6H2O": {"mw": 444.47, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 99.0}
    },
    "Y": {
        "Y(NO3)3 · 6H2O": {"mw": 383.01, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 99.0}
    },
    "Li": {
        "LiNO3 (anhydrous)": {"mw": 68.946, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 99.0},
        "LiCH3COO · 2H2O": {"mw": 102.02, "phase": "solid", "hygro": False, "type": "acetate", "default_assay": 99.0}
    },
    "Al": {
        "Al(NO3)3 · 9H2O": {"mw": 375.13, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 98.5},
        "Al(OiPr)3 (Aluminum Isopropoxide)": {"mw": 204.24, "phase": "solid", "type": "alkoxide", "default_assay": 98.0}
    },
    "P": {
        "(NH4)H2PO4": {"mw": 115.03, "phase": "solid", "hygro": False, "type": "phosphate", "default_assay": 99.0}
    },
    "Cu": {
        "Cu(NO3)2 · 3H2O": {"mw": 241.60, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 99.0}
    },
    "Zn": {
        "Zn(NO3)2 · 6H2O": {"mw": 297.49, "phase": "solid", "hygro": True, "type": "nitrate", "default_assay": 98.0}
    }
}

ORGANICS_INFO = {
    "Citric Acid (anhydrous)": {"mw": 192.124},
    "Citric Acid (monohydrate)": {"mw": 210.140},
    "Ethylene Glycol": {"mw": 62.070, "density": 1.1132}
}

# -----------------------------------------------------------------------------
# ENGINE MATHEMATICS
# -----------------------------------------------------------------------------
def parse_formula(formula_str):
    cleaned = re.sub(r'[\s_\{\}\$]', '', formula_str)
    matches = re.findall(r'([A-Z][a-z]?)([0-9]*\.?[0-9]*)', cleaned)
    comp = {}
    for elem, coeff in matches:
        if not elem:
            continue
        val = float(coeff) if coeff else 1.0
        comp[elem] = comp.get(elem, 0.0) + val
    return comp

def calculate_molecular_weight(comp):
    total = 0.0
    for el, c in comp.items():
        if el not in ATOMIC_WEIGHTS:
            raise ValueError(f"Element '{el}' is not in the atomic weight database.")
        total += c * ATOMIC_WEIGHTS[el]
    return total

def derive_pechini_chemistry(comp, app_domain):
    has_alkoxide = False
    has_carbonate_trap = False
    weighted_k1 = 0.0
    total_cations = 0.0

    for el, c in comp.items():
        if el == "O":
            continue
        total_cations += c
        spec = CATION_SPECS.get(el, {"valency": 3, "log_k1": 6.0, "alkoxide": False, "carbonate_trap": False})
        if spec.get("alkoxide", False):
            has_alkoxide = True
        if spec.get("carbonate_trap", False):
            has_carbonate_trap = True
        weighted_k1 += c * spec.get("log_k1", 6.0)

    avg_log_k1 = weighted_k1 / total_cations if total_cations > 0 else 6.0

    if has_alkoxide:
        rec_ca_ratio = 3.0 if avg_log_k1 < 10 else 3.5
        rec_eg_ratio = 4.0
        rec_ph = "7.0 – 7.5"
        mechanism = "Alkoxide solvolysis regime (requires high ligand density to suppress premature hydrolysis into TiO2/ZrO2 flocs)"
    elif app_domain == "SOFC/SOEC (High Catalytic Surface Area)":
        rec_ca_ratio = 1.5
        rec_eg_ratio = 2.5
        rec_ph = "5.5 – 6.5 (Strict: pH > 7.5 precipitates Co/Fe hydroxides!)"
        mechanism = "Fuel-lean combustion regime (minimizes carbon retention for high specific surface area)"
    elif has_carbonate_trap:
        rec_ca_ratio = 2.5
        rec_eg_ratio = 3.5
        rec_ph = "6.0 – 6.5"
        mechanism = "Alkaline-earth de-carbonation regime (requires thermal decarbonation to split BaCO3/SrCO3 intermediates)"
    else:
        rec_ca_ratio = 2.0
        rec_eg_ratio = 3.0
        rec_ph = "6.0 – 7.0"
        mechanism = "Standard transition-metal esterification"

    return rec_ca_ratio, rec_eg_ratio, rec_ph, mechanism, has_alkoxide, has_carbonate_trap

# -----------------------------------------------------------------------------
# USER INTERFACE
# -----------------------------------------------------------------------------
col1, col2 = st.columns([1.2, 1.8], gap="large")

with col1:
    st.subheader("1. Target Formulation & Sintering")
    
    app_domain = st.selectbox(
        "Application Domain:",
        [
            "SOFC/SOEC (High Catalytic Surface Area)",
            "Solid-State Battery Electrolytes / Cathodes",
            "AEM / PEM / Water Splitting Electrocatalysts"
        ]
    )

    preset = st.selectbox(
        "Validated Benchmarks or Custom:",
        [
            "LaCo0.8Fe0.2O3 (SOFC Cathode / OER)",
            "La0.6Sr0.4Co0.2Fe0.8O3 (LSCF-6428 SOFC)",
            "La0.8Sr0.2MnO3 (LSM-82 Perovskite)",
            "BaZr0.8Y0.2O3 (BZY-20 Proton Conductor)",
            "Ce0.9Gd0.1O1.95 (GDC-10 SOFC Electrolyte)",
            "Li0.40La0.60TiO3 (LLTO Battery Electrolyte)",
            "Li7La3Zr2O12 (LLZO Garnet Battery Electrolyte)",
            "NiFe2O4 (Spinel Water Splitting Catalyst)",
            "Custom Target Formulation"
        ]
    )

    preset_defaults = {
        "LaCo0.8Fe0.2O3": ("LaCo0.8Fe0.2O3", 900),
        "La0.6Sr0.4Co0.2Fe0.8O3": ("La0.6Sr0.4Co0.2Fe0.8O3", 950),
        "La0.8Sr0.2MnO3": ("La0.8Sr0.2MnO3", 1000),
        "BaZr0.8Y0.2O3": ("BaZr0.8Y0.2O3", 1250),
        "Ce0.9Gd0.1O1.95": ("Ce0.9Gd0.1O1.95", 800),
        "Li0.40La0.60TiO3": ("Li0.40La0.60TiO3", 850),
        "Li7La3Zr2O12": ("Li7La3Zr2O12", 1100),
        "NiFe2O4": ("NiFe2O4", 750)
    }

    form_default, temp_default = ("LaCo0.8Fe0.2O3", 900)
    for k, (f, t) in preset_defaults.items():
        if k in preset:
            form_default, temp_default = f, t
            break

    formula_str = st.text_input("Stoichiometric Formula Unit:", value=form_default)
    batch_mass_g = st.number_input("Target Sintered Ceramic Mass (g):", min_value=0.1, max_value=500.0, value=10.0, step=0.5)
    sinter_temp_c = st.slider("Final Sintering / Calcination Temp (°C):", 500, 1450, temp_default, 25)

    comp_prelim = parse_formula(formula_str)
    has_li = "Li" in comp_prelim
    if has_li:
        if sinter_temp_c <= 750:
            auto_li = 0.0
        elif sinter_temp_c <= 900:
            auto_li = 5.0 * ((sinter_temp_c - 750) / 150.0)
        else:
            auto_li = 5.0 + 10.0 * ((sinter_temp_c - 900) / 300.0)
        li_excess_pct = st.number_input("Volatile Li Loss Compensation (%):", 0.0, 40.0, auto_li, 1.0,
                                        help="Knudsen-Langmuir Li2O sublimation compensation factor.")
    else:
        li_excess_pct = 0.0

with col2:
    st.subheader("2. Physical Chemistry & Precursor Assay")

    try:
        comp = parse_formula(formula_str)
        mw_target = calculate_molecular_weight(comp)
        rec_ca, rec_eg, rec_ph, mechanism, has_alk, has_carb = derive_pechini_chemistry(comp, app_domain)

        m1, m2, m3 = st.columns(3)
        m1.metric("Unit MW", f"{mw_target:.3f} g/mol")
        m2.metric("Target Moles", f"{(batch_mass_g / mw_target):.5f} mol")
        m3.metric("Regime", mechanism.split()[0])
        st.info(f"💡 **Mechanistic Regime:** {mechanism} | **Safe pH:** {rec_ph}")

    except Exception as e:
        st.error(f"Invalid Formula: {e}")
        st.stop()

    st.markdown("#### Chelation & Transesterification Ratios")
    col_c1, col_c2 = st.columns(2)
    with col_c1:
        ca_type = st.selectbox("Citric Acid Form:", ["Citric Acid (anhydrous)", "Citric Acid (monohydrate)"])
        ca_ratio = st.slider("Citric Acid : Total Cation (CA:M) Molar Ratio:", 1.0, 5.0, float(rec_ca), 0.1)
    with col_c2:
        eg_ratio = st.slider("Ethylene Glycol : Citric Acid (EG:CA) Molar Ratio:", 1.5, 8.0, float(rec_eg), 0.5)

    st.markdown("#### Reagent Sourcing & Active Purity / Assay (%)")
    st.caption("Adjust Purity / Assay if your bottle has absorbed atmospheric water or is labeled < 99%.")
    
    metal_cations = [el for el in comp.keys() if el != "O"]
    selected_reagents = {}

    precursor_cols = st.columns(min(len(metal_cations), 4))
    for i, el in enumerate(metal_cations):
        col_c = precursor_cols[i % len(precursor_cols)]
        with col_c:
            if el in PRECURSOR_LIBRARY:
                opts = list(PRECURSOR_LIBRARY[el].keys())
                choice = st.selectbox(f"{el} Source:", opts, key=f"reag_{el}")
                lib_entry = PRECURSOR_LIBRARY[el][choice].copy()
                assay_val = st.number_input(f"{el} Assay / Purity (%):", min_value=50.0, max_value=105.0, 
                                            value=float(lib_entry.get("default_assay", 99.0)), step=0.5, key=f"assay_{el}")
                lib_entry["assay"] = assay_val
                lib_entry["display_name"] = choice
                selected_reagents[el] = lib_entry
            else:
                st.warning(f"Universal Fallback: {el}")
                def_mw = ATOMIC_WEIGHTS.get(el, 100.0) + (3 * 62.00) + (6 * 18.015)
                name_custom = st.text_input(f"{el} Name:", value=f"{el}(NO3)3 · 6H2O", key=f"name_{el}")
                mw_custom = st.number_input(f"{el} MW:", value=def_mw, key=f"mw_{el}")
                is_liq = st.checkbox(f"{el} Liquid?", value=False, key=f"liq_{el}")
                dens = st.number_input(f"{el} Density (g/mL):", value=1.00, key=f"den_{el}") if is_liq else 1.00
                assay_c = st.number_input(f"{el} Assay (%):", value=99.0, key=f"assay_{el}")
                selected_reagents[el] = {
                    "display_name": name_custom,
                    "mw": mw_custom,
                    "phase": "liquid" if is_liq else "solid",
                    "density": dens,
                    "assay": assay_c,
                    "type": "custom"
                }

# -----------------------------------------------------------------------------
# EXECUTION & DISPENSING
# -----------------------------------------------------------------------------
st.markdown("---")
if st.button("Generate Experimental Bench SOP", type="primary"):
    n_target = batch_mass_g / mw_target

    elem_moles = {}
    for el, coeff in comp.items():
        if el == "O":
            continue
        if el == "Li":
            elem_moles[el] = n_target * coeff * (1.0 + (li_excess_pct / 100.0))
        else:
            elem_moles[el] = n_target * coeff

    total_cation_moles = sum(elem_moles.values())

    precursor_table = []
    alkoxide_list = []

    for el, moles in elem_moles.items():
        spec = selected_reagents[el]
        purity_fraction = spec.get("assay", 100.0) / 100.0
        # Correct mass for chemical assay/purity
        reagent_mass = (moles * spec["mw"]) / purity_fraction

        if spec["phase"] == "liquid":
            vol_ml = reagent_mass / spec["density"]
            dispense_str = f"{vol_ml:.2f} mL ({reagent_mass:.3f} g)"
        else:
            dispense_str = f"{reagent_mass:.3f} g"

        if spec.get("type") == "alkoxide":
            alkoxide_list.append(f"{el} ({spec['display_name']})")

        precursor_table.append({
            "Element": el,
            "Target Reagent": spec["display_name"],
            "Theoretical Moles": f"{moles:.5f} mol",
            "Assay Applied": f"{spec.get('assay', 100.0):.1f}%",
            "Dispense Quantity": dispense_str,
            "Class": spec.get("type", "Standard").capitalize()
        })

    moles_ca = total_cation_moles * ca_ratio
    mass_ca = moles_ca * ORGANICS_INFO[ca_type]["mw"]

    moles_eg = moles_ca * eg_ratio
    mass_eg = moles_eg * ORGANICS_INFO["Ethylene Glycol"]["mw"]
    vol_eg = mass_eg / ORGANICS_INFO["Ethylene Glycol"]["density"]
    est_nh3_vol = (moles_ca * 1.2 * 17.031) / (0.90 * 0.28)

    # Physical Crucible & Foam Sizing
    min_crucible_vol_ml = max(100.0, float(batch_mass_g) * 35.0)

    organics_table = [
        {"Reagent": ca_type, "Moles": f"{moles_ca:.4f} mol", "Quantity": f"{mass_ca:.2f} g", "Role": "Chelating Complexing Agent"},
        {"Reagent": "Ethylene Glycol", "Moles": f"{moles_eg:.4f} mol", "Quantity": f"{vol_eg:.1f} mL ({mass_eg:.2f} g)", "Role": "Polyester Crosslinker"},
        {"Reagent": "Aqueous NH4OH (~28%)", "Moles": "Buffer", "Quantity": f"~{est_nh3_vol:.1f} mL (dropwise to {rec_ph})", "Role": "Carboxylic Deprotonation"}
    ]

    # ----------------- DISPLAY TABLES -----------------
    st.subheader("📋 Official Laboratory Dispensing Sheet")
    st.caption(f"Generated on {datetime.now().strftime('%Y-%m-%d %H:%M')} | Batch: {batch_mass_g:.1f}g {formula_str}")

    col_t1, col_t2 = st.columns([1.3, 1.0])
    with col_t1:
        st.markdown("#### Inorganic Metal Precursor Salts & Alkoxides")
        st.dataframe(pd.DataFrame(precursor_table), use_container_width=True, hide_index=True)
    with col_t2:
        st.markdown("#### Organic Polymerization Network")
        st.dataframe(pd.DataFrame(organics_table), use_container_width=True, hide_index=True)

    # ----------------- CRUCIBLE EXPANSION ALERT -----------------
    st.error(f"""
    🚨 **MANDATORY CRUCIBLE SIZING SAFETY REQUIREMENT:**
    * During the $200–250^\\circ\\text{{C}}$ charring step, this sol will expand **$5\\times$ to $8\\times$ in volume** as water, $\\text{{CO}}_2$, and $\\text{{NO}}_x$ boil off.
    * **Minimum Crucible Volume Required: ≥ {min_crucible_vol_ml:.0f} mL**.
    * *Warning:* Do not use standard small analytical crucibles ($30–50\\text{{ mL}}$)! Use a wide alumina saggar, boat, or tall-form $250–500\\text{{ mL}}$ beaker, or the foam will boil over into the furnace heating elements!
    """)

    # ----------------- BENCH SOP CHECKLIST -----------------
    st.subheader("🧪 Standard Operating Procedure (Step-by-Step Bench Checklist)")

    if alkoxide_list:
        solvation_instructions = f"""
- [ ] **Step 1A (Citrate-Glycol Sol):** Pour **{vol_eg/2:.1f} mL** (50%) of the Ethylene Glycol into Beaker A. Add **{mass_ca:.2f} g** of {ca_type}. Stir and heat to **50–60°C** until completely dissolved and clear.
- [ ] **Step 1B (Cooling):** Allow Beaker A to cool to **< 40°C**. *Crucial: Hot glycol accelerates violent, uncontrolled alkoxide condensation.*
- [ ] **Step 1C (Moisture-Sensitive Alkoxide Injection):** Using dry micropipette tips, slowly inject **{', '.join(alkoxide_list)}** into the vortex of Beaker A. Stir for 20 minutes. Solution must remain clear and free of precipitation[cite: 1].
- [ ] **Step 1D (Metal Salt Chelation):** In Beaker B, dissolve all non-alkoxide nitrate/acetate salts in the remaining **{vol_eg/2:.1f} mL** of Ethylene Glycol (add 15–20 mL deionized water if necessary to fully dissolve). Slowly pour Beaker B into Beaker A under constant stirring[cite: 1].
        """
    else:
        solvation_instructions = f"""
- [ ] **Step 1A (Citric Acid Dissolution):** Pour **{vol_eg/2:.1f} mL** (50%) of the Ethylene Glycol into Beaker A. Add **{mass_ca:.2f} g** of {ca_type}. Heat to **50°C** with magnetic stirring until clear.
- [ ] **Step 1B (Salt Dissolution Buffer):** In Beaker B, place all metal nitrate/acetate salts. Add the remaining **{vol_eg/2:.1f} mL** of Ethylene Glycol along with **20–30 mL of deionized water** (water buffer prevents dry nitrate clumping and hot spots). Stir at 50°C until completely dissolved.
- [ ] **Step 1C (Sol Mixing):** Pour Beaker B into Beaker A under continuous stirring at 50°C for 30 minutes to equilibrate metal-citrate coordination complexes.
        """

    st.markdown(f"""
{solvation_instructions}

- [ ] **Step 2 (pH Neutralization & Exotherm Control):**
  * Place the reaction beaker in a room-temperature or light ice-water bath.
  * Insert a calibrated pH electrode into the sol[cite: 1].
  * Slowly add aqueous $\\text{{NH}}_4\\text{{OH}}$ (~28%) dropwise while stirring vigorously[cite: 1].
  * **Target pH: {rec_ph}**.
  * *Alert:* Maintain solution temperature below **35°C** during neutralization. Exceeding $\\text{{pH}} = 8.0$ causes irreversible precipitation of transition metal hydroxides!

- [ ] **Step 3 (Polyesterification & Gelation):**
  * Move the neutralized beaker to a hotplate set to **130–140°C** with continuous magnetic stirring[cite: 1].
  * Free water and excess solvent will evaporate, driving the transesterification crosslinking reaction[cite: 1].
  * Continue heating until the sol thickens into a transparent, amber-colored, highly viscous resin (gel point)[cite: 1].

- [ ] **Step 4 (Primary Pyrolysis / Charring):**
  * Transfer the viscous gel into an alumina crucible or boat of **at least {min_crucible_vol_ml:.0f} mL capacity**[cite: 1].
  * Place in a well-ventilated drying oven or box furnace pre-heated to **220–250°C** for **1 hour**[cite: 1].
  * The resin will froth, release water vapor/$\\text{{NO}}_x$, and expand into a black, brittle solid foam[cite: 1].

- [ ] **Step 5 (Decarbonation Burnout):**
  * Transfer the charred black foam to a high-temperature box furnace under a fume hood[cite: 1].
  * Ramp at **3°C/min to 450–500°C** and dwell for **4–5 hours** in air[cite: 1].
  * This oxidizes residual hydrocarbons into $\\text{{CO}}_2$ and leaves a fine, carbon-free precursor powder[cite: 1].

- [ ] **Step 6 (Homogenization & Final Sintering):**
  * Grind the precursor powder in an agate mortar (or planetary ball mill for 30 minutes) to break agglomerates[cite: 1].
  * Sinter in an alumina crucible at **{sinter_temp_c}°C for 8–24 hours** in air[cite: 1].
  * Result: High-purity, single-phase crystalline **{formula_str}** ready for XRD and electrochemical evaluation[cite: 1].
    """)

    if has_carb:
        st.warning("⚠️ **Alkaline-Earth Carbonate Warning (Ba/Sr Present):** Citric acid pyrolysis generates $\\text{CO}_2$, which forms ultra-stable $\\text{BaCO}_3$ / $\\text{SrCO}_3$. Sintering at $\\ge 1150^\\circ\\text{C}$ with an intermediate milling step after $500^\\circ\\text{C}$ is strictly required to force complete reaction into the perovskite structure!")
