# PechiniLab — Academic Prototype

PechiniLab is a transparent, reproducible prototype for theoretical stoichiometry and precursor-quantity calculations for nominal oxide compositions.

## Important status statement

This is **not a validated synthesis protocol**. It does not predict phase purity, yield, pH, thermodynamic speciation, thermal behavior, lithium loss, or safe processing conditions. All precursor identities, molar masses, assay bases, and sources must be independently verified. Do not use generated quantities or any future process guidance without qualified review.

## Features

- Formula parsing for element symbols, decimal coefficients, and nested `()` / `[]` groups.
- Formula-mass and theoretical formula-mole calculations.
- Explicit user-defined precursor formulas and molar masses.
- Precursor mass calculations accounting for target-element atoms per precursor formula and active-compound assay.
- Liquid volume calculation when the user provides a density.
- JSON calculation records with application version, timestamp, inputs, precursor definitions, limitations, and a SHA-256 record hash.
- CSV export and unit tests.

## Not currently supported

- Hydrate-dot notation in target formulas, charges, isotope notation, phase labels, or chemical reaction equations.
- Automatic oxidation-state or charge-balance checks.
- Automatic precursor selection.
- Variable hydration / solvation models.
- Solution concentration corrections beyond the stated active-compound mass-fraction model.
- pH prediction, complexation equilibria, thermal analysis, process optimization, phase-purity or yield predictions.
- Uncertainty propagation or experimental validation.

## Install and run

Use Python 3.10 or newer.

```bash
python -m venv .venv
# macOS/Linux:
source .venv/bin/activate
# Windows PowerShell:
# .venv\Scripts\Activate.ps1

python -m pip install -r requirements.txt
streamlit run app.py
```

## Run tests

```bash
pytest -q
```

## Academic use and publication

Before describing this as research-grade or validated software:

1. Independently verify atomic weights and all precursor properties against citable sources.
2. Expand the tests and benchmark calculations against independently derived reference cases.
3. Document supported formula syntax and all assumptions.
4. Add dependency locking, continuous integration, license, authorship, citation metadata, and a versioned release.
5. For any synthesis recommendations, validate the relevant chemistry and process conditions for each material class.
6. If claiming improved materials outcomes, collect and analyze experimental evidence.

The included tests are initial software tests, not evidence that the application has been scientifically validated.
