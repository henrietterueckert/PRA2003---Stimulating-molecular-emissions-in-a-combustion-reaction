# PRA2003 - Simulating molecular emissions in a combustion reaction
Simulating molecular emissions in a combustion reaction. Analysing 5 million simulated combustion events for 12 molecular species with their counterparts.
Calculating average abundance and checking for asymmetries between paired species

**Student**

Henriette Rückert - i6397668 

## Dependencies

- **Python**: 3.8 or later

- **Dependencies**: import math, csv, statistics but no pip install needed

- **Input files**: output-Set1.txt through output-Set10.txt must be in the same directory as the script. Each line is either a header (eventNumber nParticles) or a particle row (px py pz moleculeID)

- **Run from repository**:
`week2deliverable.py`
`week3deliverable.py`
`week4deliverable.py`

- **Runtime**: several minutes.

## Data ##

**Week 2:** The input is output-Set0.txt

**Week 3:** The input is any of the output files

**Week 4:** The input is 10 text files

output-Set1.txt --> output-Set10.txt, each containing 500,000 events (5,000,000 events total)

Each event in a file has:

- A Header line: eventID  and number of molecules rows
  
- The molecule rows contain: px py pz  ad the moleculeID. These are the 3D momentum components with integer ID identifying the molecule

An empty event is a real event where nothing happened, script excludes it from the event count N

## Output files ##

**Week 2**

The output is printed directly onto console

**Week 3** 

The output is printed directly onto console

**Week 4**
The output is printed directly onto console and into csv files:

`subsample_results.csv` = Per-code results for each of the 10 files separately 

`results.csv` = Pooled average per code across all 10 files, with sub-sampling uncertainty 

`significance.csv` = Pairwise difference, correlation, and significance test for each molecule/counterpart pair 

## Method ##

**Week 2**

Reads one event from a single file, computes each particle's momentum magnitude, and prints everything to the console

**Week 3**

It only handles one file and one chosen ID per run. Excludes empty events from N, and computes the uncertainty as the standard error of the mean

**Week 4**

1. Read one file, tally every known code per (non-empty) event. Keeps a running total per code:
   
- total_count: the sum of counts, used to get the mean
- Poisson: sqrt(N) / n_events
   
2. Combine the 10 files into one final average per code and takes the uncertainty as the standard deviation of the 10 per-file averages (sub-sampling method)
   
3. Pairwise test: for each code pair, the script computes the difference and asymmetry per file first, then the standard deviation of those per-file differences is the uncertainty on the difference. correlation_r is reported alongside as a diagnostic.

## How does the sub-sampling method work? ##

The 5 million event sample is split into 10 sub-samples of 500,000 events each and a central value is taken from the average over the full pooled sample. The statistical uncertainty is the standard deviation of the 10 sub-sample averages. This is used instead of propagating each molecule's uncertainty independently, since we don't know correlation between molecules measured in the same events.

# Questions
# What are the average counts of each molecular species and their statistical uncertainties?

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

# Is there any asymmetry between the normal and the variant molecule?

For each molecule/counterpart pair, the difference and asymmetry are calculated per sub-sample first, then the standard deviation across sub-samples gives the uncertainty on the difference (n_sigma). This accounts for correlation between the pair. Pairs with n_sigma > 3 are significant. z_sem and z_std are shown for comparison  as they  would combine each side's uncertainty independently.

| ID pair | Molecule (vs counterpart) | Difference | Uncertainty | n_sigma | Asymmetry [%] | Significant? |
| --- | --- | --- | --- | --- | --- | --- | 
| 211 vs -211 | Carbon monoxide (vs Carbon-13 monoxide) | 0.0323 | ± 0.005 |  7.15 | 0.081 | Yes |
| 321 vs -321 | Nitric oxide (vs Ionised NO) | 0.0057 | ± 0.003 | 1.73 | 0.114 | No |
| 2212 vs -2212 | Water (vs Heavy water) | 0.0239 | ± 0.002 | 10.06 | 0.998 | Yes |
| 3122 vs -3122 | Methane (vs Methyl ion) | 0.0049 | ± 0.001| 8.41 |  0.894 | Yes |
| 3312 vs -3312 | Ethylene (vs Ionised ethylene) | 0.0004 | ± 0.000 | 0.90 | 0.562 | No |
| 3334 vs -3334 | Ozone (vs Superoxide anion) | 0.0000 | ± 0.000 | 0.51 |  1.496 | No |

3 of 6 pairs show a significant asymmetry (> 3 sigma): Carbon monoxide, Water, and Methane.
  
# Is there any asymmetry as a function of their momentum?
## Plots

**`asymmetry_vs_P.png`** has one panel per pair, showing:
- Blue points: the asymmetry in each bin. The vertical bar is the std of the sub-samples. The horizontal bar is the bin width
- Orange band: the overall asymmetry ± its std
- Green slope: the straight-line fit, with the slope given in the legend
- Dashed grey line: A = 0

**`spectrum_dNdP.png`** has one panel per pair, with the normal molecule in blue and the variant in red, on a log scale. The open last bin (P > 8) is not included, because it has no real width.


## Notes and limitations

- Rare pairs: some pairs, such as ozone / superoxide, have very few molecules. Their asymmetries have large uncertainties.
- Bins not plotted: a bin is left off the asymmetry plot if fewer than 2 sub-samples have molecules in it, or if its std is 0.
- What a zero slope means: a slope consistent with zero does not prove the asymmetry is independent of P (it could rise and then fall). The per-bin comparison with the overall asymmetry looks at that
