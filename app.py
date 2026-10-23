import streamlit as st
import pandas as pd
from pechini_core import calculate_llto_pechini, deconvolute_icp_ppm

st.set_page_config(page_title="Pechini Sol-Gel Toolbox", page_icon="⚗️", layout="wide")

st.title("⚗️ Pechini Sol-Gel & Perovskite Synthesis Toolbox")
st.markdown("Automating laboratory stoichiometry, alkoxide volume dosing, and ICP-MS site-occupancy analysis for ceramic solid electrolytes.")

tab1, tab2 = st.tabs(["🧪 1. Lab Precursor & Chelation Calculator", "🔬 2. ICP-MS Elemental Deconvolution"])

with tab1:
    st.subheader("Pechini Synthesis Formulation: Li_x La_y TiO_3")
    c1, c2, c3 = st.columns(3)
    with c1:
        x_li = st.number_input("Li target (x)", 0.01, 1.0, 0.40, 0.01)
        y_la = st.number_input("La target (y)", 0.01, 1.0, 0.60, 0.01)
        target_mass = st.number_input("Target Dry Mass (g)", 0.5, 500.0, 10.0, 0.5)
    with c2:
        la_hydrated = st.checkbox("La(NO3)3 is Hexahydrate (6H2O)", value=True)
        ca_monohydrate = st.checkbox("Citric Acid is Monohydrate (H2O)", value=False)
        st.caption("Ti precursor: Ti(IV) butoxide liquid (density = 1.00 g/mL)")
    with c3:
        ca_ratio = st.slider("Citric Acid : Metal Cation Ratio", 1.0, 5.0, 3.0, 0.5)
        eg_factor = st.slider("Citric Acid : Ethylene Glycol Ratio (1:N)", 1, 10, 8, 1)

    if st.button("Generate Lab Recipe", type="primary"):
        res = calculate_llto_pechini(
            x_li=x_li, y_la=y_la, target_mass_g=target_mass,
            la_hydrated=la_hydrated, citric_acid_monohydrate=ca_monohydrate,
            ca_to_metal_ratio=ca_ratio, ca_to_eg_molar_ratio=(1.0 / eg_factor)
        )
        st.success(f"Recipe for **{res['formula']}** | Molecular Weight: **{res['mw']:.2f} g/mol**")
        
        col_a, col_b = st.columns(2)
        with col_a:
            st.markdown("### 🧱 Inorganic Precursors")
            p_data = [
                {"Precursor": "LiNO3", "Amount to Weigh": f"{res['mass_lino3']:.3f} g", "Moles": f"{res['moles_li']:.5f}"},
                {"Precursor": "La(NO3)3 · 6H2O" if la_hydrated else "La(NO3)3", "Amount to Weigh": f"{res['mass_lano3']:.3f} g", "Moles": f"{res['moles_la']:.5f}"},
                {"Precursor": "Ti(OBu)4 (Liquid)", "Amount to Pipette": f"{res['vol_tibuo4']:.2f} mL ({res['mass_tibuo4']:.3f} g)", "Moles": f"{res['moles_ti']:.5f}"},
            ]
            st.dataframe(pd.DataFrame(p_data), use_container_width=True, hide_index=True)

        with col_b:
            st.markdown("### 🧬 Organic Complexing Agents")
            o_data = [
                {"Component": "Citric Acid", "Amount to Weigh": f"{res['mass_ca']:.2f} g", "Moles": f"{res['moles_ca']:.4f}"},
                {"Component": "Ethylene Glycol", "Amount to Measure": f"{res['vol_eg']:.1f} mL ({res['mass_eg']:.2f} g)", "Moles": f"{res['moles_eg']:.4f}"},
                {"Component": "Ammonia (NH3)", "Amount": "Add dropwise until pH = 7.0", "Role": "pH buffer & chelation stabilizer"}
            ]
            st.dataframe(pd.DataFrame(o_data), use_container_width=True, hide_index=True)

with tab2:
    st.subheader("ICP-MS Cation Occupancy Deconvolution")
    st.markdown("Input elemental PPM readings to calculate true stoichiometry normalized to Ti = 1.0.")
    preset = st.radio("Load Sample Data:", ["Synthesized LLTO (LLTO-CHR-1)", "Commercial Benchmark (LLTO-COM-1)", "Custom"], horizontal=True)
    if preset == "Synthesized LLTO (LLTO-CHR-1)":
        init_la, init_li, init_ti = 418300.0, 16877.0, 263200.0
    elif preset == "Commercial Benchmark (LLTO-COM-1)":
        init_la, init_li, init_ti = 410400.0, 15643.0, 272800.0
    else:
        init_la, init_li, init_ti = 415000.0, 16000.0, 265000.0

    k1, k2, k3 = st.columns(3)
    la_val = k1.number_input("La (ppm)", value=init_la, step=1000.0)
    li_val = k2.number_input("Li (ppm)", value=init_li, step=500.0)
    ti_val = k3.number_input("Ti (ppm)", value=init_ti, step=1000.0)

    if st.button("Analyze ICP Data", type="primary"):
        df_res = deconvolute_icp_ppm({"La": la_val, "Li": li_val, "Ti": ti_val})
        st.dataframe(df_res.style.format({
            "PPM (mg/kg)": "{:,.0f}",
            "Wt (%)": "{:.2f}%",
            "Moles / 100g": "{:.5f}",
            "Mol (%)": "{:.2f}%",
            "Stoichiometry (Ti=1.0)": "{:.2f}"
        }), use_container_width=True)
        
        li_calc = df_res.loc[df_res["Element"] == "Li", "Stoichiometry (Ti=1.0)"].values[0]
        la_calc = df_res.loc[df_res["Element"] == "La", "Stoichiometry (Ti=1.0)"].values[0]
        st.metric("Deconvoluted Formula", f"Li_{{{li_calc:.2f}}} La_{{{la_calc:.2f}}} TiO_3", delta=f"Target: Li_0.40 La_0.60 TiO_3")
