import streamlit as st
import pandas as pd
import numpy as np
import re

# -----------------------------------------------------------------------------
# PAGE CONFIGURATION
# -----------------------------------------------------------------------------
st.set_page_config(
    page_title="PechiniLab | Advanced Electroceramics Suite",
    page_icon="⚗️",
    layout="wide"
)

st.title("⚗️ PechiniLab: Universal Sol-Gel Suite for Electroceramics")
st.markdown(
    "A thermodynamics- and coordination-chemistry-grounded engine for **SOFC/SOEC Cathodes & Electrolytes**, "
    "**Solid-State Battery Conductors**, and **OER/AEM/PEM Electrocatalysts**."
)

# -----------------------------------------------------------------------------
# THERMODYNAMIC & CHEMICAL KNOWLEDGE BASE
# -----------------------------------------------------------------------------
ATOMIC_WEIGHTS = {
    "H": 1.008, "Li": 6.941, "B": 10.811, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "Si": 28.085, "P": 30.974,
    "S": 32.065, "K": 39.098, "Ca": 40.078, "Sc": 44.956, "Ti": 47.867, "V": 50.942,
    "Cr": 51.996, "Mn": 54.938, "Fe": 55.845, "Co": 58.933, "Ni": 58.693, "Cu": 63.546,
    "Zn": 65.380, "Ga": 69.723, "Ge": 72.630, "Sr": 87.620, "Y": 88.906, "Zr": 91.224,
    "Nb": 92.906, "Mo": 95.950, "Ru": 101.07, "Rh": 102.91, "Pd": 106.42, "Ag": 107.87,
    "In": 114.82, "Sn": 118.71, "Sb": 121.76, "Ba": 137.33, "La": 138.91, "Ce": 140.12,
    "Pr": 140.91, "Nd": 144.24, "Sm": 150.36, "Eu": 151.96, "Gd": 157.25, "Tb": 158.93,
    "Dy": 162.50, "Ho": 164.93, "Er": 167.26, "Tm": 168.93, "Yb": 173.05, "Lu": 174.97,
    "Hf": 178.49, "Ta": 180.95, "W": 183.84, "Ir": 192.22, "Pt": 195.08, "Pb": 207.20,
    "Bi": 208.98
}

# Empirical chelation stability log K1 with citrate and default valency
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

# Precursor specifications with precise hydration and densities
PRECURSOR_LIBRARY = {
    "La": {
        "La(NO3)3 · 6H2O": {"mw": 433.01, "phase": "solid", "hygro": True, "type": "nitrate"},
        "La(CH3COO)3 · 1.5H2O": {"mw": 343.04, "phase": "solid", "hygro": False, "type": "acetate"}
    },
    "Co": {
        "Co(NO3)2 · 6H2O": {"mw": 291.03, "phase": "solid", "hygro": True, "type": "nitrate"},
        "Co(CH3COO)2 · 4H2O": {"mw": 249.08, "phase": "solid", "hygro": False, "type": "acetate"}
    },
    "Fe": {
        "Fe(NO3)3 · 9H2O": {"mw": 404.00, "phase": "solid", "hygro": True, "type": "nitrate"},
        "Fe(CH3COO)2": {"mw": 173.93, "phase": "solid", "hygro": False, "type": "acetate"},
        "FeCl3 (anhydrous)": {"mw": 162.20, "phase": "solid", "hygro": True, "type": "chloride"}
    },
    "Sr": {
        "Sr(NO3)2 (anhydrous)": {"mw": 211.63, "phase": "solid", "hygro": False, "type": "nitrate"},
        "Sr(CH3COO)2": {"mw": 205.71, "phase": "solid", "hygro": False, "type": "acetate"}
    },
    "Ba": {
        "Ba(NO3)2 (anhydrous)": {"mw": 261.34, "phase": "solid", "hygro": False, "type": "nitrate"},
        "Ba(CH3COO)2": {"mw": 255.42, "phase": "solid", "hygro": False, "type": "acetate"}
    },
    "Ni": {
        "Ni(NO3)2 · 6H2O": {"mw": 290.79, "phase": "solid", "hygro": True, "type": "nitrate"},
        "Ni(CH3COO)2 · 4H2O": {"mw": 248.86, "phase": "solid", "hygro": False, "type": "acetate"}
    },
    "Mn": {
        "Mn(NO3)2 · 4H2O": {"mw": 251.01, "phase": "solid", "hygro": True, "type": "nitrate"},
        "Mn(CH3COO)2 · 4H2O": {"mw": 245.09, "phase": "solid", "hygro": False, "type": "acetate"}
    },
    "Ti": {
        "Ti(OBu)4 (Titanium(IV) Butoxide, 97%)": {"mw": 340.35, "phase": "liquid", "density": 1.000, "type": "alkoxide"},
        "Ti(OiPr)4 (Titanium(IV) Isopropoxide, 97%)": {"mw": 284.22, "phase": "liquid", "density": 0.960, "type": "alkoxide"}
    },
    "Zr": {
        "Zr(OPr)4 (Zirconium Propoxide 70% in PrOH)": {"mw": 327.57 / 0.70, "phase": "liquid", "density": 1.044, "type": "alkoxide"},
        "ZrO(NO3)2 · xH2O (Zirconyl Nitrate)": {"mw": 231.23, "phase": "solid", "hygro": True, "type": "nitrate"}
    },
    "Ce": {
        "Ce(NO3)3 · 6H2O": {"mw": 434.22, "phase": "solid", "hygro": True, "type": "nitrate"}
    },
    "Gd": {
        "Gd(NO3)3 · 6H2O": {"mw": 451.36, "phase": "solid", "hygro": True, "type": "nitrate"}
    },
    "Sm": {
        "Sm(NO3)3 · 6H2O": {"mw": 444.47, "phase": "solid", "hygro": True, "type": "nitrate"}
    },
    "Y": {
        "Y(NO3)3 · 6H2O": {"mw": 383.01, "phase": "solid", "hygro": True, "type": "nitrate"}
    },
    "Li": {
        "LiNO3 (anhydrous)": {"mw": 68.946, "phase": "solid", "hygro": True, "type": "nitrate"},
        "LiCH3COO · 2H2O": {"mw": 102.02, "phase": "solid", "hygro": False, "type": "acetate"}
    },
    "Al": {
        "Al(NO3)3 · 9H2O": {"mw": 375.13, "phase": "solid", "hygro": True, "type": "nitrate"},
        "Al(OiPr)3 (Aluminum Isopropoxide)": {"mw": 204.24, "phase": "solid", "type": "alkoxide"}
    },
    "P": {
        "(NH4)H2PO4": {"mw": 115.03, "phase": "solid", "hygro": False, "type": "phosphate"}
    }
}

ORGANICS_INFO = {
    "Citric Acid (anhydrous)": {"mw": 192.124},
    "Citric Acid (monohydrate)": {"mw": 210.140},
    "Ethylene Glycol": {"mw": 62.070, "density": 1.1132}
}

# -----------------------------------------------------------------------------
# CHEMICAL & FORMULA PARSER ENGINE
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
    """
    Literature-grounded determination of CA:M, EG:CA, and pH based on cation speciation.
    """
    has_alkoxide = False
    has_carbonate_trap = False
    max_k1 = 0.0
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
        k1 = spec.get("log_k1", 6.0)
        max_k1 = max(max_k1, k1)
        weighted_k1 += c * k1

    avg_log_k1 = weighted_k1 / total_cations if total_cations > 0 else 6.0

    # Literature Rules for Chelation & Polymerization Ratios:
    if has_alkoxide:
        # Ti4+/Zr4+ alkoxides require high ligand density to suppress hydrolysis into TiO2/ZrO2 flocs
        rec_ca_ratio = 3.0 if avg_log_k1 < 10 else 3.5
        rec_eg_ratio = 4.0
        rec_ph = "7.0 – 7.5"
        mechanism = "Alkoxide-chelated polymeric network (hydrolysis-suppression regime)"
    elif app_domain == "SOFC/SOEC (High Catalytic Surface Area)":
        # Pure nitrates (e.g. La-Sr-Co-Fe-O) require lean CA to prevent excessive carbon retention
        rec_ca_ratio = 1.5
        rec_eg_ratio = 2.5
        rec_ph = "5.5 – 6.5 (Prevent precipitation of Fe(OH)3 / Co(OH)2)"
        mechanism = "Fuel-lean combustion / porous resin regime"
    elif has_carbonate_trap:
        # Ba/Sr-rich systems form persistent BaCO3; requires high CA excess and high-T burnout
        rec_ca_ratio = 2.5
        rec_eg_ratio = 3.5
        rec_ph = "6.0 – 6.5"
        mechanism = "Alkaline-earth de-carbonation regime"
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
    st.subheader("1. Material System & Domain")
    
    app_domain = st.selectbox(
        "Application Domain:",
        [
            "SOFC/SOEC (High Catalytic Surface Area)",
            "Solid-State Battery Electrolytes / Cathodes",
            "AEM / PEM / Water Splitting Electrocatalysts"
        ]
    )

    preset = st.selectbox(
        "Literature Benchmarks or Custom:",
        [
            "LaCo0.8Fe0.2O3 (SOFC Cathode / OER)",
            "La0.6Sr0.4Co0.2Fe0.8O3 (LSCF-6428 SOFC)",
            "La0.8Sr0.2MnO3 (LSM-82 Perovskite)",
            "BaZr0.8Y0.2O3 (BZY-20 Proton Conductor)",
            "Ce0.9Gd0.1O1.95 (GDC-10 SOFC Electrolyte)",
            "Li0.40La0.60TiO3 (LLTO Battery Electrolyte)",
            "Li7La3Zr2O12 (LLZO Garnet Battery Electrolyte)",
            "NiFe2O4 (Spinel Anion Exchange Catalyst)",
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

    formula_str = st.text_input("Stoichiometric Formula:", value=form_default)
    batch_mass_g = st.number_input("Target Sintered Ceramic Mass (g):", min_value=0.1, max_value=500.0, value=10.0, step=0.5)
    sinter_temp_c = st.slider("Calcination / Sintering Temperature (°C):", 500, 1400, temp_default, 25)

    # Volatilization physics
    comp_prelim = parse_formula(formula_str)
    has_li = "Li" in comp_prelim
    if has_li:
        if sinter_temp_c <= 750:
            auto_li = 0.0
        elif sinter_temp_c <= 900:
            auto_li = 5.0 * ((sinter_temp_c - 750) / 150.0)
        else:
            auto_li = 5.0 + 10.0 * ((sinter_temp_c - 900) / 300.0)
        li_excess_pct = st.number_input("Volatile Li Loss Compensation (%):", 0.0, 40.0, auto_li, 1.0)
    else:
        li_excess_pct = 0.0

with col2:
    st.subheader("2. Physical Chemistry & Sol-Gel Parameters")

    try:
        comp = parse_formula(formula_str)
        mw_target = calculate_molecular_weight(comp)
        rec_ca, rec_eg, rec_ph, mechanism, has_alk, has_carb = derive_pechini_chemistry(comp, app_domain)

        m1, m2, m3 = st.columns(3)
        m1.metric("Molecular Weight", f"{mw_target:.3f} g/mol")
        m2.metric("Target Moles", f"{(batch_mass_g / mw_target):.5f} mol")
        m3.metric("Regime", mechanism.split()[0])
        st.info(f"💡 **Chemical Mechanism:** {mechanism} | Recommended pH: **{rec_ph}**")

    except Exception as e:
        st.error(f"Invalid Formula: {e}")
        st.stop()

    st.markdown("#### Chelation & Transesterification Ratios")
    col_c1, col_c2 = st.columns(2)
    with col_c1:
        ca_type = st.selectbox("Citric Acid Precursor:", ["Citric Acid (anhydrous)", "Citric Acid (monohydrate)"])
        ca_ratio = st.slider(
            "Citric Acid : Total Cation (CA:M) Molar Ratio:",
            min_value=1.0, max_value=5.0, value=float(rec_ca), step=0.1,
            help="Higher values guarantee chelation stability but require longer burnout."
        )
    with col_c2:
        eg_ratio = st.slider(
            "Ethylene Glycol : Citric Acid (EG:CA) Molar Ratio:",
            min_value=1.5, max_value=8.0, value=float(rec_eg), step=0.5,
            help="Controls crosslinking density vs. combustion exotherm."
        )

    st.markdown("#### Reagent Sourcing (Automated Fallback for Any Cation)")
    metal_cations = [el for el in comp.keys() if el != "O"]
    selected_reagents = {}

    precursor_cols = st.columns(min(len(metal_cations), 4))
    for i, el in enumerate(metal_cations):
        col_c = precursor_cols[i % len(precursor_cols)]
        with col_c:
            if el in PRECURSOR_LIBRARY:
                opts = list(PRECURSOR_LIBRARY[el].keys())
                choice = st.selectbox(f"{el} Precursor:", opts, key=f"reag_{el}")
                selected_reagents[el] = PRECURSOR_LIBRARY[el][choice]
                selected_reagents[el]["display_name"] = choice
            else:
                st.warning(f"Fallback for {el}")
                def_mw = ATOMIC_WEIGHTS.get(el, 100.0) + (3 * 62.00) + (6 * 18.015)
                name_custom = st.text_input(f"{el} Name:", value=f"{el}(NO3)3 · 6H2O", key=f"name_{el}")
                mw_custom = st.number_input(f"{el} MW:", value=def_mw, key=f"mw_{el}")
                is_liq = st.checkbox(f"{el} Liquid?", value=False, key=f"liq_{el}")
                dens = st.number_input(f"{el} Density (g/mL):", value=1.00, key=f"den_{el}") if is_liq else 1.00
                selected_reagents[el] = {
                    "display_name": name_custom,
                    "mw": mw_custom,
                    "phase": "liquid" if is_liq else "solid",
                    "density": dens,
                    "type": "custom"
                }

# -----------------------------------------------------------------------------
# EXECUTION & BENCH SOP
# -----------------------------------------------------------------------------
st.markdown("---")
if st.button("Generate Scientifically Validated Protocol", type="primary"):
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
        reagent_mass = moles * spec["mw"]

        if spec["phase"] == "liquid":
            vol_ml = reagent_mass / spec["density"]
            dispense_str = f"**{vol_ml:.2f} mL** micropipette ({reagent_mass:.3f} g)"
        else:
            dispense_str = f"**{reagent_mass:.3f} g** balance"

        if spec.get("type") == "alkoxide":
            alkoxide_list.append(f"{el} ({spec['display_name']})")

        precursor_table.append({
            "Element": el,
            "Reagent": spec["display_name"],
            "Cation Moles": f"{moles:.5f} mol",
            "Dispense Quantity": dispense_str,
            "Classification": spec.get("type", "Standard").capitalize()
        })

    moles_ca = total_cation_moles * ca_ratio
    mass_ca = moles_ca * ORGANICS_INFO[ca_type]["mw"]

    moles_eg = moles_ca * eg_ratio
    mass_eg = moles_eg * ORGANICS_INFO["Ethylene Glycol"]["mw"]
    vol_eg = mass_eg / ORGANICS_INFO["Ethylene Glycol"]["density"]

    organics_table = [
        {"Reagent": ca_type, "Moles": f"{moles_ca:.4f} mol", "Quantity": f"**{mass_ca:.2f} g**", "Role": "Chelating Complexing Agent"},
        {"Reagent": "Ethylene Glycol", "Moles": f"{moles_eg:.4f} mol", "Quantity": f"**{vol_eg:.1f} mL** ({mass_eg:.2f} g)", "Role": "Polyester Crosslinker"},
        {"Reagent": "Aqueous NH4OH (~28%)", "Moles": "Buffer", "Quantity": f"Add dropwise to {rec_ph}", "Role": "Deprotonation Stabilizer"}
    ]

    st.subheader("📋 Experimental Dispensing Matrix")
    c_t1, c_t2 = st.columns([1.3, 1.0])
    with c_t1:
        st.markdown("#### Inorganic Precursor Salts & Alkoxides")
        st.dataframe(pd.DataFrame(precursor_table), use_container_width=True, hide_index=True)
    with c_t2:
        st.markdown("#### Polymeric Complexation Network")
        st.dataframe(pd.DataFrame(organics_table), use_container_width=True, hide_index=True)

    # -------------------------------------------------------------------------
    # SCIENTIFIC BENCH SOP
    # -------------------------------------------------------------------------
    st.subheader("🧪 Literature-Grounded Standard Operating Procedure (SOP)")

    if alkoxide_list:
        seq_text = f"""
        * **Sequential Alkoxide Chelation Protocol (Crucial):**
          1. Dissolve **{mass_ca:.2f} g** of {ca_type} into **50% of the EG ({vol_eg/2:.1f} mL)** at 60°C until completely clear.
          2. Cool this acidic glycol solution to **< 45°C**.
          3. Under dry Ar or a vigorous vortex, slowly inject the alkoxides: **{', '.join(alkoxide_list)}**. The alkoxide must exchange ligands with citrate/glycol *before* encountering any water of hydration.
          4. Separately dissolve the hydrated salts in the remaining EG, then merge both solutions slowly.
        """
    else:
        seq_text = f"""
        * **Direct Solvation Protocol (All Nitrates / Acetates):**
          1. Dissolve **{mass_ca:.2f} g** of {ca_type} in **50% of the EG ({vol_eg/2:.1f} mL)** at 50°C.
          2. Dissolve all metal salts into the remaining EG with ~20 mL of deionized water to accelerate dissolution.
          3. Combine both beakers under continuous stirring at 50°C.
        """

    st.markdown(f"""
    1. **Pre-Dissolution & Chelation Architecture:**
       {seq_text}

    2. **Deprotonation & pH Neutralization Window:**
       * Slowly add aqueous $\\text{{NH}}_4\\text{{OH}}$ dropwise while monitoring with a calibrated pH probe.
       * **Target pH: {rec_ph}**.
       * *Scientific Warning:* Exceeding $\\text{{pH}} = 8.0$ causes catastrophic precipitation of divalent transition metals (e.g., $\\text{{Co(OH)}}_2$, $\\text{{Fe(OH)}}_3$), ruining cation stoichiometry!

    3. **Condensation & Resin Gelation:**
       * Heat the neutral solution to **130–140°C** with constant magnetic stirring.
       * Water will evaporate, driving the transesterification reaction between the carboxyl groups of citrate and the hydroxyl groups of ethylene glycol.
       * Maintain heating until the sol transitions into an amber, transparent, highly viscous polyester resin (gel point).

    4. **Primary Pyrolysis / Charring:**
       * Transfer the gel into an alumina crucible.
       * Heat to **220–250°C** for 1 hour in a drying oven or box furnace. The resin will expand into a voluminous, brittle black polymer foam.

    5. **Calcination & High-Temperature Sintering:**
       * **Organic Burnout:** Ramp at $3^\\circ\\text{{C/min}}$ to **$450–500^\\circ\\text{{C}}$ for 4–5 hours** in air.
       * **Milling:** Pulverize the charred powder in an agate mortar or planetary ball mill to eliminate hard agglomerates.
       * **Final Sintering:** Heat in a high-temperature box furnace to **{sinter_temp_c}°C for 8–24 hours** to complete crystallization of pure **{formula_str}**.
    """)

    if has_carb:
        st.warning("⚠️ **Carbonate Trapping Alert (Barium/Strontium System):** Alkaline-earth cations readily react with decomposing citrate fuels to form ultra-stable $\\text{BaCO}_3$ / $\\text{SrCO}_3$. Sintering at $>1150^\\circ\\text{C}$ or performing an intermediate milling step is strictly required to force solid-state reaction with B-site oxides!")
