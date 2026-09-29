# PRA2003 - Stimulating molecular emissions in a combustion reaction

**Student**

Henriette Rückert - i6397668 
## Dependencies

- **Python**: 3.8 or later

- **Dependencies**: import math, csv, statistics but no pip install needed

- **Input files**: output-Set1.txt through output-Set10.txt must be in the same directory as the script. Each line is either a header (eventNumber nParticles) or a particle row (px py pz moleculeID)

- **Run**: full analysis.py

- **Runtime**: a few minutes total across all 10 files, depending on machine

## Questions ##
**What are the average counts of each molecular species and their statistical uncertainties?**

| ID | Molecule | Average per event | Statistical uncertainty |
| --- | --- | --- | --- |
| 211 | Carbon monoxide | 19.949508 | ± 0.032746 |
| -211 | Carbon-13 monoxide | 19.917211 | ± 0.031873 |
| 321 | Nitric oxide | 2.509148 | ± 0.004774 |
| -321 | Ionised NO | 2.503457 | ± 0.005504 |
| 2212 | Water | 1.208034 | ± 0.001896 |
| -2212 | Heavy water (D2O) | 1.184161 | ± 0.002412 |
| 3122 | Methane (CH4) | 0.276599 | ± 0.001074 |
| -3122 | Methyl ion (CH3-) | 0.271696 | ± 0.000985 |
| 3312 | Ethylene (C2H4) | 0.039441 | ± 0.000284 |
| -3312 | Ionised ethylene (C2H3-) | 0.039000 | ± 0.000402 |
| 3334 | Ozone (O3) | 0.001187 | ± 0.000042 |
| -3334 | Superoxide anion (O2-) | 0.001152 | ± 0.000051 |

**Is there any asymmetry between the normal and the variant molecule?**

| ID pair | Molecule (vs counterpart) | Difference | Uncertainty | Correlation r | n_sigma | z_sem | z_std | Asymmetry [%] | Significant? |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 211 vs -211 | Carbon monoxide (vs Carbon-13 monoxide) | 0.0323 | ± 0.005 | 0.991 | 7.15 | 2.23 | 0.71 | 0.081 | Yes |
| 321 vs -321 | Nitric oxide (vs Ionised NO) | 0.0057 | ± 0.003 | 0.804 | 1.73 | 2.47 | 0.78 | 0.114 | No |
| 2212 vs -2212 | Water (vs Heavy water) | 0.0239 | ± 0.002 | 0.414 | 10.06 | 24.60 | 7.78 | 0.998 | Yes |
| 3122 vs -3122 | Methane (vs Methyl ion) | 0.0049 | ± 0.001 | 0.843 | 8.41 | 10.64 | 3.36 | 0.894 | Yes |
| 3312 vs -3312 | Ethylene (vs Ionised ethylene) | 0.0004 | ± 0.000 | 0.019 | 0.90 | 2.83 | 0.90 | 0.562 | No |
| 3334 vs -3334 | Ozone (vs Superoxide anion) | 0.0000 | ± 0.000 | -0.090 | 0.51 | 1.68 | 0.53 | 1.496 | No |

3 of 6 pairs show a significant asymmetry (> 3 sigma): Carbon monoxide, Water, and Methane.

Carbon monoxide has the strongest correlation of any pair (r = 0.991). Because of that,  z-scores that treat the two sides as independent (z_sem = 2.23, z_std = 0.71) stay below the significance threshold but n_sigma (7.15), which looks at how the difference itself varies across sub-samples instead of combining each side's uncertainty separately, shows asymmetry.
  
**Is there any asymmetry as a function of their momentum?**

To be answered next week
