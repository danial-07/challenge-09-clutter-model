# Data Challenge 09 – Growth and Yield with the Clutter Model

Nine permanent plots (C531, 43 measurements, ages 2.58–6.37 years) at four site classes: S = 20 m (2 plots), 22 m (3), 24 m (3) and 26 m (**1 plot**). All numbers below come from `scripts/01_clutter_model.R`.

## 1. Fitted equations

**Equation 1 – basal-area projection** (82 pairs of measurements)

ln(B2) = (I1/I2)·ln(B1) + 2.8109·(1 − I1/I2) + 0.02330·(1 − I1/I2)·S

| Coefficient | Estimate | Std. error | t |
|---|---|---|---|
| α0 | 2.8109 | 0.1914 | 14.7 |
| α1 | 0.02330 | 0.00832 | 2.8 |

R² (on ln B2) = 0.929; Syx = 0.0486 in ln units = **5.1 %** of mean B2. (On the transformed variable the R² is 0.990; the lower value above is the honest one for predicting B2.)

**Equation 2 – static volume** (43 measurements)

ln(V) = 1.3050 − 1.1319·(1/I) + 0.01602·S + 1.2551·ln(B)

| Coefficient | Estimate | Std. error | t |
|---|---|---|---|
| β0 | 1.3050 | 0.1782 | 7.3 |
| β1 | −1.1319 | 0.2641 | −4.3 |
| β2 | 0.01602 | 0.00586 | 2.7 |
| β3 | 1.2551 | 0.0730 | 17.2 |

R² = 0.991; Syx = 0.0491 in ln units = **5.0 %** of mean V. All signs are biologically sensible: more basal area, better site and older age give more volume.

**Residuals** (`outputs/residuals_equations_1_and_2.png`):

![Residuals of Equations 1 and 2](../outputs/residuals_equations_1_and_2.png)

With only 43–82 points the residuals show no strong trend. The 82 pairs share the same measurements (a plot with 5 measurements gives 10 pairs), so they are **not independent**: the standard errors and R² of Equation 1 are somewhat optimistic.

## 2. Compatible volume equation (Equation 3)

Equation 3 was used in its **sequential form** (B2 from Equation 1, then V2 from Equation 2 at age I2), which the assignment accepts. The script also writes it analytically and confirms both forms give identical values:

ln(V2) = β0 + β1/I2 + β2·S + β3·[(I1/I2)·ln(B1) + (α0 + α1·S)(1 − I1/I2)]

**Note on the assignment PDF.** The printed Equation 3 (with β0 and β2·S multiplied by 1 − I1/I2, and without ln B1) is not what results from substituting Equation 1 into Equation 2. Used literally, it underestimates V2 by about 66 % on our data, so we used the derived form. This should be confirmed against Chapter 11.

**Does the projection reproduce the observed volumes?** Over the 82 pairs, the mean bias is −0.3 % and the RMSE is **7.6 %** of mean V2. Individual errors range from −12 % to +25 %. Leave-one-plot-out (refitting both equations without the plot being predicted) gives an RMSE of 9.6 %, a fairer indication of how the model behaves on a new plot. By site (in-sample) the RMSE is 4.3 % (S=20), 6.9 % (S=22), 9.4 % (S=24) and 4.5 % (S=26), with V2 underestimated by about 5 % at S=22 and overestimated by about 6 % at S=24, so the S effect is smooth while the plots are not perfectly ordered by S. In the single-plot test (plot 2-2, S=22) the projection from 2.62 to 3.64 years is almost exact (86.99 vs 86.96 m³/ha), but projections to 6.30 years underestimate the observed 171.6 m³/ha by 7–10 %.

![Projected vs observed V2](../outputs/eq3_observed_vs_projected.png)

## 3. Yield prognosis

**Initial condition:** age 2.58 years (youngest observed) and basal area 7.79 m²/ha (mean of the first measurement of the nine plots), the same for the four sites, so differences come only from S. Age step 0.5 year, up to 10 years; ages above 6.37 years are extrapolation.

![Yield curves](../outputs/yield_curve.png)

| Age (years) | S = 20 | S = 22 | S = 24 | S = 26 |
|---|---|---|---|---|
| 3.58 | 74.8 | 78.5 | 82.4 | 86.5 |
| 5.58 | 124.6 | 132.8 | 141.5 | 150.7 |
| 9.58 (extrapolated) | 182.4 | 196.6 | 211.9 | 228.3 |

(Volume, m³/ha, from `yield_prognosis.csv`.)

**Differences between sites grow with age.** The gap between S = 26 and S = 20 is 4 m³/ha at 2.6 years, 26 m³/ha at 5.6 years and 46 m³/ha at 9.6 years; the ratio V(26)/V(20) rises from 1.10 to 1.21 to 1.25. This follows from the model: the site index enters Equation 1 as an asymptote (α1·S), so better sites keep accumulating basal area for longer.

## 4. ICA, IMA and technical rotation age

![ICA and IMA by site](../outputs/ica_ima_chart.png)

| Site (m) | Rotation age (ICA = IMA) | IMA at rotation (m³/ha/yr) | Volume at rotation (m³/ha) |
|---|---|---|---|
| 20 | 5.63 | 22.3 | 125.6 |
| 22 | 5.78 | 23.7 | 137.1 |
| 24 | 5.93 | 25.3 | 149.7 |
| 26 | 6.08 | 26.9 | 163.5 |

In every site ICA starts well above IMA, falls steadily, and crosses the almost flat IMA curve; after that ICA stays below IMA. The crossing moves **later by about 0.15 year (≈ 2 months) per 2 m of site index**, while the maximum IMA rises from 22 to 27 m³/ha/yr. All four crossings fall inside the observed age range (≤ 6.37 years), but close to its upper end.

## 5. Limitations (please read before using these numbers)

- **The ordering of the rotation ages is not robust.** It depends on the initial basal area. With a common B0 (above), better sites rotate later. If each site starts from its own plots' mean B0 (6.5, 6.9, 8.9 and 9.8 m²/ha), the rotation ages become 6.22, 6.19, 5.49 and 5.35 years: the order reverses (see `initial_condition_sensitivity.csv`). The four values are within about 0.9 year of each other under either choice, so we do **not** conclude that the technical rotation differs between sites.
- **Sample size.** 9 plots; the S = 26 m class rests on a single plot and S = 20 m on two. The S coefficients (t ≈ 2.7–2.8) are driven by very few independent units.
- **Age range.** Data cover 2.6–6.4 years. The rotation ages (5.6–6.1 years) are at the edge of the data, and everything beyond 6.4 years (dashed lines in the graphs) is extrapolation of a log-linear model. The early peak of IMA (around 6 years) is earlier than commonly reported for commercial eucalyptus, which suggests the model's saturation behaviour should not be trusted at older ages.
- **In-sample fit.** The equations were fitted and evaluated on the same plots; the leave-one-plot-out RMSE (9.6 %) is the more realistic figure.
- **Back-transformation bias.** Predictions come from log models without a correction factor; the factor exp(Syx²/2) is ≈ 1.001 here, so the effect is negligible.
