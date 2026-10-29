"""
Unit tests for PechiniLab core stoichiometry and chemical parser.
"""

import pytest
from app import parse_formula, molecular_weight, calculate_pechini_synthesis, PrecursorConfig

def test_hydrate_parsing():
    comp = parse_formula("La(NO3)3*6H2O")
    assert comp["La"] == 1.0
    assert comp["N"] == 3.0
    assert comp["H"] == 12.0
    assert comp["O"] == 15.0  # 9 + 6

def test_dot_notation():
    comp = parse_formula("Al(NO3)3·9H2O")
    assert comp["Al"] == 1.0
    assert comp["N"] == 3.0
    assert comp["H"] == 18.0
    assert comp["O"] == 18.0

def test_complex_perovskite():
    comp = parse_formula("Li0.33La0.55TiO3")
    assert pytest.approx(comp["Li"], 0.001) == 0.33
    assert pytest.approx(comp["La"], 0.001) == 0.55
    assert pytest.approx(comp["Ti"], 0.001) == 1.0
    assert pytest.approx(comp["O"], 0.001) == 3.0

def test_garnet_llzo():
    comp = parse_formula("Li7La3Zr2O12")
    assert comp["Li"] == 7.0
    assert comp["La"] == 3.0
    assert comp["Zr"] == 2.0
    assert comp["O"] == 12.0

def test_batch_calculation():
    precursors = {
        "Li": PrecursorConfig("Li", "LiNO3", "LiNO3", 68.946, 1.0, 99.0),
        "La": PrecursorConfig("La", "La(NO3)3*6H2O", "La(NO3)3*6H2O", 433.01, 1.0, 99.0),
        "Ti": PrecursorConfig("Ti", "TTIP", "Ti(OCH(CH3)2)4", 284.22, 1.0, 99.0, phase="liquid", density_g_ml=0.96)
    }
    calc = calculate_pechini_synthesis(
        "Li0.33La0.55TiO3",
        10.0,
        precursors,
        {"Li": 10.0},
        ca_to_metal_ratio=2.0,
        eg_to_ca_ratio=4.0
    )
    assert calc["target_mass_g"] == 10.0
    assert calc["total_metal_moles"] > 0
    assert calc["organics"]["ca_mass_g"] > 0
    assert calc["organics"]["eg_volume_ml"] > 0
