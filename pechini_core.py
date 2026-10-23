"""
Pechini Sol-Gel Synthesis Precursor & Chelation Calculator + ICP-MS Deconvolution.
Based on experimental workflows for LLTO perovskites (Aristotle University of Thessaloniki).
"""
import pandas as pd

MW_DATA = {
    "Li": 6.941,
    "La": 138.905,
    "Ti": 47.867,
    "O": 15.999,
    "LiNO3": 68.946,
    "La(NO3)3_6H2O": 433.01,
    "La(NO3)3_anhydrous": 324.91,
    "Ti(OBu)4": 340.35,          # Titanium(IV) butoxide
    "CitricAcid_anhydrous": 192.124,
    "CitricAcid_monohydrate": 210.14,
    "EthyleneGlycol": 62.07,
}

DENSITY_DATA = {
    "Ti(OBu)4": 1.00,            # g/mL
    "EthyleneGlycol": 1.1132     # g/mL
}

def calculate_llto_pechini(
    x_li=0.40,
    y_la=0.60,
    target_mass_g=10.0,
    la_hydrated=True,
    citric_acid_monohydrate=False,
    ca_to_metal_ratio=3.0,
    ca_to_eg_molar_ratio=1.0 / 8.0,
):
    mw_llto = (x_li * MW_DATA["Li"]) + (y_la * MW_DATA["La"]) + (1.0 * MW_DATA["Ti"]) + (3.0 * MW_DATA["O"])
    moles_target = target_mass_g / mw_llto

    moles_li = moles_target * x_li
    moles_la = moles_target * y_la
    moles_ti = moles_target * 1.0
    total_cations = moles_li + moles_la + moles_ti

    # Inorganic precursors
    mass_lino3 = moles_li * MW_DATA["LiNO3"]
    la_mw = MW_DATA["La(NO3)3_6H2O"] if la_hydrated else MW_DATA["La(NO3)3_anhydrous"]
    mass_lano3 = moles_la * la_mw

    mass_tibuo4 = moles_ti * MW_DATA["Ti(OBu)4"]
    vol_tibuo4 = mass_tibuo4 / DENSITY_DATA["Ti(OBu)4"]

    # Organics (Citric Acid & Ethylene Glycol)
    moles_ca = total_cations * ca_to_metal_ratio
    ca_mw = MW_DATA["CitricAcid_monohydrate"] if citric_acid_monohydrate else MW_DATA["CitricAcid_anhydrous"]
    mass_ca = moles_ca * ca_mw

    moles_eg = moles_ca / ca_to_eg_molar_ratio
    mass_eg = moles_eg * MW_DATA["EthyleneGlycol"]
    vol_eg = mass_eg / DENSITY_DATA["EthyleneGlycol"]

    return {
        "formula": f"Li_{{{x_li:.2f}}}La_{{{y_la:.2f}}}TiO_3",
        "mw": mw_llto,
        "mass_lino3": mass_lino3,
        "moles_li": moles_li,
        "mass_lano3": mass_lano3,
        "moles_la": moles_la,
        "mass_tibuo4": mass_tibuo4,
        "vol_tibuo4": vol_tibuo4,
        "moles_ti": moles_ti,
        "mass_ca": mass_ca,
        "moles_ca": moles_ca,
        "mass_eg": mass_eg,
        "vol_eg": vol_eg,
        "moles_eg": moles_eg,
    }

def deconvolute_icp_ppm(ppm_dict):
    data = []
    for element in ["La", "Li", "Ti"]:
        ppm = ppm_dict.get(element, 0.0)
        wt_percent = ppm / 10000.0  # 1% = 10,000 ppm
        mol_100g = wt_percent / MW_DATA[element]
        data.append({
            "Element": element,
            "PPM (mg/kg)": ppm,
            "Wt (%)": wt_percent,
            "Moles / 100g": mol_100g
        })

    df = pd.DataFrame(data)
    total_moles = df["Moles / 100g"].sum()
    df["Mol (%)"] = (df["Moles / 100g"] / total_moles) * 100.0

    # Normalize to Ti = 1.0 (B-site cation)
    moles_ti = df.loc[df["Element"] == "Ti", "Moles / 100g"].values[0]
    df["Stoichiometry (Ti=1.0)"] = df["Moles / 100g"] / moles_ti
    return df
