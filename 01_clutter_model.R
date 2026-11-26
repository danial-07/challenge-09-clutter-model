# =============================================================================
# 01_clutter_model.R
# Data Challenge 09 - Growth and yield with the Clutter (1963) model
#
# HOW TO RUN
#   1. Install once:  install.packages("readxl")
#   2. setwd() to the folder "challenge-09-clutter-model" (the one that has the
#      folders scripts/, data/, outputs/).
#   3. source("scripts/01_clutter_model.R")
#
# WHAT IT DOES (the 8 steps of the assignment)
#   1. read data, build pairs of measurements          -> outputs/pairs_base.csv
#   2. fit Equation 1 (basal-area projection)          -> outputs/clutter_equations.csv
#   3. fit Equation 2 (static volume)                  -> outputs/clutter_equations.csv
#   4. test the compatible volume equation (Eq. 3)     -> outputs/eq3_test_*.csv
#   5. simulate yield for the 4 sites                  -> outputs/yield_prognosis.csv
#   6. ICA and IMA, technical rotation age             -> outputs/technical_rotation.csv
#   7. graphs                                          -> outputs/*.png
# =============================================================================

# ----------------------------- CONFIG (edit here) ----------------------------
DATA_FILE      <- "data/Clutter_data_original.xlsx"      # input file
OUTPUT_DIR     <- "outputs"                              # output folder
FUNCTIONS_FILE <- "scripts/00_functions_challenge09.R"   # functions file

AGE_STEP           <- 0.5   # age step of the simulation, in years (must divide 1)
MAX_SIMULATION_AGE <- 10    # last simulated age. Ages above the oldest observed
                            # age are flagged "extrapolated" in the output.
INITIAL_CONDITION  <- "common"  # "common"   = same B0 for the 4 sites (mean B of the
                                #              first measurement of all plots)
                                # "per_site" = B0 = mean first-measurement B of
                                #              the plots of each site
TEST_PLOT          <- "2-2"     # plot used for the single-plot test of Equation 3
# -----------------------------------------------------------------------------

source(FUNCTIONS_FILE)
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
out <- function(name) file.path(OUTPUT_DIR, name)

# ---- Step 1: data and pairs --------------------------------------------------
data <- read_clutter(DATA_FILE)
pairs <- build_pairs(data)
cat("Measurements:", nrow(data), "| plots:", length(unique(data$plot_id)),
    "| pairs:", nrow(pairs), "\n")
cat("Plots per site:\n"); print(tapply(data$plot_id, data$S,
                                       function(x) length(unique(x))))
write.csv(pairs, out("pairs_base.csv"), row.names = FALSE)

# ---- Steps 2 and 3: fit Equations 1 and 2 ------------------------------------
eq1 <- fit_eq1(pairs)
eq2 <- fit_eq2(data)
alpha <- eq1$alpha     # c(alpha0, alpha1)
beta <- eq2$beta     # c(beta0, beta1, beta2, beta3)

eq_table <- rbind(
  data.frame(equation = "Eq1_basal_area", group = "coefficient",
             name = c("alpha0", "alpha1"), value = alpha, std_error = eq1$std_error),
  data.frame(equation = "Eq1_basal_area", group = "metric",
             name = names(eq1$metrics), value = unname(eq1$metrics), std_error = NA),
  data.frame(equation = "Eq2_static_volume", group = "coefficient",
             name = c("beta0", "beta1", "beta2", "beta3"), value = beta,
             std_error = eq2$std_error),
  data.frame(equation = "Eq2_static_volume", group = "metric",
             name = names(eq2$metrics), value = unname(eq2$metrics), std_error = NA)
)
write.csv(eq_table, out("clutter_equations.csv"), row.names = FALSE)
cat("\nEquation 1:\n"); print(round(eq1$metrics, 4)); print(round(alpha, 5))
cat("\nEquation 2:\n"); print(round(eq2$metrics, 4)); print(round(beta, 5))

# ---- Step 4: Equation 3 (compatible volume) ----------------------------------
# 4a. the sequential form and the analytical form must be identical
pairs$V2_proj <- project_V(pairs$I1, pairs$B1, pairs$I2, pairs$S, alpha, beta)
V2_analytic  <- project_V_analytic(pairs$I1, pairs$B1, pairs$I2, pairs$S, alpha, beta)
stopifnot(isTRUE(all.equal(pairs$V2_proj, V2_analytic)))
cat("\nEq. 3: sequential and analytical forms give identical results. OK\n")

pairs$B2_proj <- project_B(pairs$I1, pairs$B1, pairs$I2, pairs$S, alpha)
pairs$rel_error_V2_pct <- 100 * (pairs$V2_proj - pairs$V2) / pairs$V2
error_summary <- function(obs, pred) {
  c(bias_pct = 100 * mean(obs - pred) / mean(obs),
    RMSE_pct = 100 * sqrt(mean((obs - pred)^2)) / mean(obs))
}
cat("Eq. 3 on all pairs (bias %, RMSE %):\n")
print(round(error_summary(pairs$V2, pairs$V2_proj), 2))

# 4b. the form printed in the PDF, for comparison (diagnostic)
V2_pdf <- project_V_pdf_literal(pairs$I1, pairs$B1, pairs$I2, pairs$S, alpha, beta)
cat("Eq. 3 exactly as printed in the PDF (bias %, RMSE %):\n")
print(round(error_summary(pairs$V2, V2_pdf), 2))

# 4c. single-plot test requested in the assignment
cat("\nSingle-plot test, plot", TEST_PLOT, ":\n")
test_plot <- pairs[pairs$plot == TEST_PLOT,
                       c("plot", "S", "I1", "I2", "B1", "B2", "B2_proj",
                         "V2", "V2_proj", "rel_error_V2_pct")]
print(round(test_plot[, -1], 2))
write.csv(pairs, out("eq3_test_pairs.csv"), row.names = FALSE)
write.csv(test_plot, out("eq3_test_plot.csv"), row.names = FALSE)

# ---- Step 5: initial condition and simulation --------------------------------
I0 <- min(data$age)                      # youngest observed age
MAX_OBS_AGE <- max(data$age)           # oldest observed age
first_meas <- do.call(rbind, lapply(split(data, data$plot_id),
                                  function(g) g[which.min(g$age), ]))
B0_common <- mean(first_meas$B)
B0_site <- tapply(first_meas$B, first_meas$S, mean)
sites   <- sort(unique(data$S))
B0_for <- function(S, mode) if (mode == "common") B0_common else unname(B0_site[as.character(S)])

cat("\nInitial age I0 =", round(I0, 3), "| common B0 =", round(B0_common, 3), "m2/ha\n")

prog <- do.call(rbind, lapply(sites, function(S)
  simulate_prognosis(S, I0, B0_for(S, INITIAL_CONDITION), MAX_SIMULATION_AGE, AGE_STEP,
                   alpha, beta, MAX_OBS_AGE)))
rownames(prog) <- NULL

# sanity check: repeated one-step projection == one direct projection from I0
V_direct <- mapply(function(S, I) project_V(I0, B0_for(S, INITIAL_CONDITION), I, S, alpha, beta),
                   prog$site, prog$age)
stopifnot(isTRUE(all.equal(prog$volume_m3ha, V_direct)))
write.csv(prog, out("yield_prognosis.csv"), row.names = FALSE)

# ---- Step 6: technical rotation age (ICA = IMA) ------------------------------
rot <- do.call(rbind, lapply(sites, function(S) {
  B0 <- B0_for(S, INITIAL_CONDITION)
  id <- technical_rotation_age(S, I0, B0, alpha, beta)
  V  <- project_V(I0, B0, id, S, alpha, beta)
  data.frame(site = S, technical_rotation_age = id, volume_at_rotation = V,
             ima_at_rotation = V / id)
}))
write.csv(rot, out("technical_rotation.csv"), row.names = FALSE)
cat("\nTechnical rotation age by site:\n"); print(round(rot, 2))

# sensitivity: how much does the rotation age depend on the initial B0?
sens <- do.call(rbind, lapply(sites, function(S) {
  data.frame(site = S,
             B0_common = B0_common,
             rotation_B0_common = technical_rotation_age(S, I0, B0_common, alpha, beta),
             B0_per_site = unname(B0_site[as.character(S)]),
             rotation_B0_per_site = technical_rotation_age(S, I0, unname(B0_site[as.character(S)]),
                                                          alpha, beta))
}))
write.csv(sens, out("initial_condition_sensitivity.csv"), row.names = FALSE)
cat("\nSensitivity to the initial condition:\n"); print(round(sens, 2))

# ---- Step 7: graphs (base R, PNG) -------------------------------------------
colors <- c("#86b6ef", "#3987e5", "#1c5cab", "#0d366b")   # blue ramp: light = low S
names(colors) <- as.character(sites)

# 7a. yield curves, all sites overlaid
png(out("yield_curve.png"), width = 1200, height = 750, res = 150)
par(mar = c(4.2, 4.5, 3, 5))
plot(NA, xlim = c(I0, MAX_SIMULATION_AGE + 1), ylim = c(0, max(prog$volume_m3ha) * 1.05),
     xlab = "Age (years)", ylab = "Volume (m3/ha)",
     main = "Yield curve by site (dashed = extrapolation)")
grid(nx = NA, ny = NULL, col = "#e1e0d9", lty = 1)
abline(v = MAX_OBS_AGE, col = "#898781")
for (S in sites) {
  g  <- prog[prog$site == S, ]
  g1 <- g[g$age <= MAX_OBS_AGE, ]
  g2 <- g[g$age >= max(g1$age), ]
  lines(g1$age, g1$volume_m3ha, col = colors[as.character(S)], lwd = 3)
  lines(g2$age, g2$volume_m3ha, col = colors[as.character(S)], lwd = 3, lty = 2)
  text(max(g$age) + 0.1, tail(g$volume_m3ha, 1), paste0("S=", S), pos = 4, xpd = TRUE, cex = 0.8)
}
legend("topleft", legend = paste0("S = ", sites, " m"), col = colors, lwd = 3, bty = "n")
dev.off()

# 7b. ICA and IMA, one panel per site, technical rotation age marked
png(out("ica_ima_chart.png"), width = 1500, height = 1050, res = 150)
par(mfrow = c(2, 2), mar = c(4, 4.2, 2.5, 1))
yl <- c(0, max(prog$ica_m3ha, prog$ima_m3ha, na.rm = TRUE) * 1.05)
for (S in sites) {
  g  <- prog[prog$site == S, ]
  ir <- rot$technical_rotation_age[rot$site == S]
  plot(NA, xlim = c(I0, MAX_SIMULATION_AGE), ylim = yl, xlab = "Age (years)",
       ylab = "m3/ha/year", main = paste0("Site S = ", S, " m"))
  rect(MAX_OBS_AGE, yl[1], MAX_SIMULATION_AGE + 1, yl[2], col = "#f0efec", border = NA)
  grid(nx = NA, ny = NULL, col = "#e1e0d9", lty = 1)
  lines(g$age, g$ica_m3ha, col = "#2a78d6", lwd = 2.5)
  lines(g$age, g$ima_m3ha, col = "#eb6834", lwd = 2.5)
  abline(v = ir, lty = 3)
  text(ir + 0.1, yl[2] * 0.93, paste0("technical rotation\n", round(ir, 2), " years"),
       pos = 4, cex = 0.8)
  if (S == sites[1]) legend("topright", legend = c("ICA", "IMA"),
                             col = c("#2a78d6", "#eb6834"), lwd = 2.5, bty = "n")
}
dev.off()

# 7c. residuals of Equations 1 and 2
png(out("residuals_equations_1_and_2.png"), width = 1500, height = 600, res = 150)
par(mfrow = c(1, 2), mar = c(4.2, 4.5, 3, 1))
plot(eq1$ln_fitted, log(pairs$B2) - eq1$ln_fitted, pch = 19, col = "#2a78d6",
     xlab = "Fitted value", ylab = "Residual (log scale)", main = "Equation 1 (ln B2)")
abline(h = 0, col = "#898781")
plot(eq2$ln_fitted, log(data$V) - eq2$ln_fitted, pch = 19, col = "#2a78d6",
     xlab = "Fitted value", ylab = "Residual (log scale)", main = "Equation 2 (ln V)")
abline(h = 0, col = "#898781")
dev.off()

# 7d. Equation 3: projected vs observed V2
png(out("eq3_observed_vs_projected.png"), width = 825, height = 750, res = 150)
par(mar = c(4.2, 4.5, 3, 1))
lim <- max(pairs$V2, pairs$V2_proj) * 1.05
plot(pairs$V2, pairs$V2_proj, pch = 19, col = "#2a78d6", xlim = c(0, lim), ylim = c(0, lim),
     xlab = "Observed V2 (m3/ha)", ylab = "V2 projected by Eq. 3 (m3/ha)",
     main = "Equation 3: projected vs observed")
abline(0, 1, col = "#898781")
dev.off()

cat("\nDone. Files written to:", normalizePath(OUTPUT_DIR), "\n")
