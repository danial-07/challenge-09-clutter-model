# Data Challenge 09 – Growth and Yield with the Clutter Model

Forest Mensuration course. Fits the Clutter (1963) system of equations to nine permanent plots of project C531 (eucalyptus), simulates the yield of the four site classes (S = 20, 22, 24, 26 m), and finds the technical rotation age (where ICA = IMA) of each site.

## Folder contents

| Path | What it is |
|---|---|
| `data/Clutter_data_original.xlsx` | Input: 43 measurements, 9 plots, ages 2.58–6.37 years |
| `scripts/00_functions_challenge09.R` | All functions (reading, pairs, the two regressions, the three equations, simulation, ICA/IMA, rotation age) |
| `scripts/01_clutter_model.R` | Main script: CONFIG block at the top, runs the 8 steps of the assignment |
| `outputs/clutter_equations.csv` | Coefficients (with standard errors) and fit metrics of Equations 1 and 2 |
| `outputs/yield_prognosis.csv` | `site, age, basal_area_m2ha, volume_m3ha, ica_m3ha, ima_m3ha` (+ `extrapolated`) |
| `outputs/technical_rotation.csv` | Technical rotation age per site |
| `outputs/initial_condition_sensitivity.csv` | How the rotation age changes with the initial basal area |
| `outputs/eq3_test_pairs.csv`, `eq3_test_plot.csv` | Projected vs observed V2 (all pairs, and one plot) |
| `outputs/yield_curve.png` | Yield curves, one per site |
| `outputs/ica_ima_chart.png` | ICA and IMA per site, rotation age marked |
| `outputs/residuals_equations_1_and_2.png`, `eq3_observed_vs_projected.png` | Residual analysis and Equation 3 check |
| `report/report_challenge09.md` | Short report with the graphs and the discussion |

## How to run

```r
install.packages("readxl")        # once
setwd("path/to/challenge-09-clutter-model")
source("scripts/01_clutter_model.R")
```

All editable settings are in the CONFIG block at the top of `01_clutter_model.R`.

## Method in plain words

1. **Pairs.** For each plot, every pair of ages (I1 < I2) gives one row with B1, B2 and V2. 4 or 5 measurements per plot give 82 pairs.
2. **Equation 1 – basal-area projection:** `ln(B2) = (I1/I2)·ln(B1) + α0·(1 − I1/I2) + α1·(1 − I1/I2)·S`. The coefficient of ln(B1) is fixed by the model at I1/I2, so the term is moved to the left side and a regression through the origin is fitted: `ln(B2) − (I1/I2)·ln(B1) = α0·x1 + α1·x2`.
3. **Equation 2 – static volume:** `ln(V) = β0 + β1/I + β2·S + β3·ln(B)`, fitted on all 43 measurements.
4. **Equation 3 – compatible volume:** computed **sequentially** (B2 from Equation 1, then V2 from Equation 2 at age I2). The script also writes it analytically and checks that both forms give the same numbers.
5. **Simulation:** Equation 1 is applied repeatedly from the initial condition; Equation 2 gives V at each age.
6. **ICA(I) = V(I) − V(I−1)**, **IMA(I) = V(I)/I.** The rotation age is where ICA crosses IMA, found exactly with `uniroot()`.

## Initial condition and age step (as requested)

- **Initial age:** 2.58 years (the youngest age observed).
- **Initial basal area:** 7.79 m²/ha, the mean B of the first measurement of the nine plots (their ages are 2.58–2.69 years). The **same B0 is used for all four sites**, so the differences between sites come only from S. A second option (`INITIAL_CONDITION <- "per_site"`: each site starts from its own plots' mean B) is available, and its effect is shown in `initial_condition_sensitivity.csv`.
- **Age step:** 0.5 year. The step must divide 1 year so that ICA = V(I) − V(I−1) can be computed.
- **Horizon:** up to 10 years. Rows older than 6.37 years (oldest observed age) are marked `extrapolated = TRUE` and drawn dashed/shaded in the graphs.

## Important notes

- All file, folder, function, variable and column names in this package are in English. The input workbook was renamed to `Clutter_data_original.xlsx`; its content is unchanged, so the column headers **inside** the workbook are still the original Portuguese ones. The script renames them by position, so this does not matter.

- **The Equation 3 printed in the assignment PDF is not the result of substituting Equation 1 into Equation 2.** The PDF writes `ln(V2) = β0(1 − I1/I2) + β1/I2 + β2·S·(1 − I1/I2) + β3·ln(B2)`. Substituting correctly gives `ln(V2) = β0 + β1/I2 + β2·S + β3·[(I1/I2)·ln(B1) + (α0 + α1·S)(1 − I1/I2)]`. The assignment allows the sequential implementation, which is what we used. The literal PDF form underestimates V2 by about 66 % on the data (the script prints this). Please confirm the notation in Chapter 11 of the textbook.
- The PNG graphs in `outputs/` are provided with the package; running the R script overwrites them with the R versions (same file names).
- Sample limits: 9 plots, only **one plot at S = 26 m**, ages 2.6–6.4 years; the 82 pairs are not independent.
