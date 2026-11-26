# =============================================================================
# 00_functions_challenge09.R
# Data Challenge 09 - Growth and yield with the Clutter (1963) model
#
# This file holds ONLY functions (no data, no plots). It is loaded by
# 01_clutter_model.R with source(). Every function has a short comment that
# says what goes in, what comes out, and why it exists.
#
# Notation (same as the assignment PDF):
#   I1, I2  = initial and final age (years)         B1, B2 = basal area (m2/ha)
#   V, V2   = volume (m3/ha)                        S      = site index (m)
#   alpha0, alpha1 = coefficients of Equation 1       beta0..beta3 = Equation 2
# =============================================================================


# -----------------------------------------------------------------------------
# 1. READING AND ORGANISING THE DATA
# -----------------------------------------------------------------------------

# Read Clutter_data_original.xlsx.
# The first 8 columns are renamed BY POSITION (the original header has accents,
# e.g. "Talhao" with a tilde, which can break on some computers).
# Output: data.frame with columns project, stand, plot, age, B, Hd, V, S
#         + plot_id (stand-plot, e.g. "2-3"), sorted by plot and age.
read_clutter <- function(path) {
  d <- as.data.frame(readxl::read_excel(path))
  names(d)[1:8] <- c("project", "stand", "plot", "age", "B", "Hd", "V", "S")
  d$plot_id <- paste(d$stand, d$plot, sep = "-")
  d <- d[order(d$plot_id, d$age), ]
  rownames(d) <- NULL
  d
}

# Build ALL pairs of measurements (I1 < I2) inside each plot.
# A plot measured 4 times gives 6 pairs, 5 times gives 10 pairs.
# (Same idea as the algebraic-difference pairs of Challenge 07.)
# Output columns: plot, S, I1, I2, B1, B2, V2
build_pairs <- function(data) {
  rows_list <- lapply(split(data, data$plot_id), function(g) {
    g <- g[order(g$age), ]
    if (nrow(g) < 2) return(NULL)
    ij <- utils::combn(nrow(g), 2)          # row 1 = first index, row 2 = second
    data.frame(plot = g$plot_id[ij[1, ]],
               S  = g$S[ij[1, ]],
               I1 = g$age[ij[1, ]], I2 = g$age[ij[2, ]],
               B1 = g$B[ij[1, ]],     B2 = g$B[ij[2, ]],
               V2 = g$V[ij[2, ]],
               stringsAsFactors = FALSE)
  })
  pairs <- do.call(rbind, rows_list)
  rownames(pairs) <- NULL
  pairs
}


# -----------------------------------------------------------------------------
# 2. FITTING THE TWO REGRESSIONS
# -----------------------------------------------------------------------------

# Equation 1 (basal-area projection):
#   ln(B2) = (I1/I2) ln(B1) + alpha0 (1 - I1/I2) + alpha1 (1 - I1/I2) S
# The coefficient of ln(B1) is NOT free: it is fixed at I1/I2 by the model.
# So we move that term to the left side and fit a regression THROUGH THE ORIGIN:
#   y = ln(B2) - (I1/I2) ln(B1)  =  alpha0 * x1 + alpha1 * x2
#   x1 = (1 - I1/I2)   x2 = (1 - I1/I2) * S
# Output: list(alpha, std_error, metrics, ln_fitted, model)
fit_eq1 <- function(pairs) {
  r  <- pairs$I1 / pairs$I2
  dd <- data.frame(y  = log(pairs$B2) - r * log(pairs$B1),
                   x1 = 1 - r,
                   x2 = (1 - r) * pairs$S)
  m  <- lm(y ~ 0 + x1 + x2, data = dd)
  n  <- nrow(dd)
  p  <- 2
  ln_obs  <- log(pairs$B2)
  ln_fit <- r * log(pairs$B1) + unname(fitted(m))   # fitted ln(B2)
  sse     <- sum((ln_obs - ln_fit)^2)
  metrics <- c(
    n_pairs = n,
    # R2 on the ln(B2) scale (the variable we really want to predict)
    R2_lnB2 = 1 - sse / sum((ln_obs - mean(ln_obs))^2),
    # R2 on the transformed variable y (usually higher; shown for transparency)
    R2_transformed_variable = 1 - sum(resid(m)^2) / sum((dd$y - mean(dd$y))^2),
    # standard error of the estimate in ln units
    Syx_ln = sqrt(sse / (n - p)),
    # standard error of the estimate in % of mean B2, on the original scale
    Syx_percent_B2 = 100 * sqrt(sum((pairs$B2 - exp(ln_fit))^2) / (n - p)) /
                        mean(pairs$B2)
  )
  list(alpha = unname(coef(m)),
       std_error = unname(summary(m)$coefficients[, "Std. Error"]),
       metrics = metrics, ln_fitted = ln_fit, model = m)
}

# Equation 2 (static volume):  ln(V) = b0 + b1 (1/I) + b2 S + b3 ln(B)
# Uses ALL measurements (it does not need two ages).
# Output: list(beta, std_error, metrics, ln_fitted, model)
fit_eq2 <- function(data) {
  dd <- data.frame(lnV = log(data$V), inv_age = 1 / data$age,
                   S = data$S, lnB = log(data$B))
  m  <- lm(lnV ~ inv_age + S + lnB, data = dd)
  n  <- nrow(dd)
  p  <- 4
  ln_fit <- unname(fitted(m))
  metrics <- c(
    n_measurements = n,
    R2_lnV = 1 - sum(resid(m)^2) / sum((dd$lnV - mean(dd$lnV))^2),
    Syx_ln = sqrt(sum(resid(m)^2) / (n - p)),
    Syx_percent_V = 100 * sqrt(sum((data$V - exp(ln_fit))^2) / (n - p)) /
                       mean(data$V)
  )
  list(beta = unname(coef(m)),
       std_error = unname(summary(m)$coefficients[, "Std. Error"]),
       metrics = metrics, ln_fitted = ln_fit, model = m)
}


# -----------------------------------------------------------------------------
# 3. THE THREE EQUATIONS AS FUNCTIONS (vectorised: accept vectors)
# -----------------------------------------------------------------------------

# Equation 1: project basal area from age I1 to age I2.
project_B <- function(I1, B1, I2, S, alpha) {
  r <- I1 / I2
  exp(r * log(B1) + (alpha[1] + alpha[2] * S) * (1 - r))
}

# Equation 2: volume from age, site and basal area (static equation).
static_volume <- function(I, S, B, beta) {
  exp(beta[1] + beta[2] / I + beta[3] * S + beta[4] * log(B))
}

# Equation 3 (compatible volume), SEQUENTIAL form:
# first B2 from Equation 1, then V2 from Equation 2 at age I2.
project_V <- function(I1, B1, I2, S, alpha, beta) {
  static_volume(I2, S, project_B(I1, B1, I2, S, alpha), beta)
}

# Equation 3, ANALYTICAL form (Equation 1 substituted inside Equation 2):
#   ln(V2) = b0 + b1/I2 + b2 S + b3 [ (I1/I2) ln(B1) + (a0 + a1 S)(1 - I1/I2) ]
# It must give exactly the same number as project_V(); the main script checks.
project_V_analytic <- function(I1, B1, I2, S, alpha, beta) {
  r <- I1 / I2
  exp(beta[1] + beta[2] / I2 + beta[3] * S +
      beta[4] * (r * log(B1) + (alpha[1] + alpha[2] * S) * (1 - r)))
}

# Equation 3 EXACTLY AS PRINTED IN THE ASSIGNMENT PDF (diagnostic only):
#   ln(V2) = b0 (1 - I1/I2) + b1/I2 + b2 S (1 - I1/I2) + b3 ln(B2)
# It is NOT what you get by substituting Eq. 1 into Eq. 2 (b0 and b2*S should
# appear in full, not multiplied by (1 - I1/I2)). Kept here only so the main
# script can show how far this printed form is from the observed volumes.
project_V_pdf_literal <- function(I1, B1, I2, S, alpha, beta) {
  r  <- I1 / I2
  B2 <- project_B(I1, B1, I2, S, alpha)
  exp(beta[1] * (1 - r) + beta[2] / I2 + beta[3] * S * (1 - r) + beta[4] * log(B2))
}


# -----------------------------------------------------------------------------
# 4. YIELD PROGNOSIS, ICA, IMA AND TECHNICAL ROTATION AGE
# -----------------------------------------------------------------------------

# Simulate B and V for ONE site, starting from (I0, B0).
# Equation 1 is applied repeatedly, one age step at a time (step = step).
# Rows with age > max_obs_age are flagged extrapolated = TRUE.
# Output columns: site, age, basal_area_m2ha, volume_m3ha,
#                 ica_m3ha, ima_m3ha, extrapolated
simulate_prognosis <- function(S, I0, B0, max_age, step, alpha, beta, max_obs_age) {
  ages <- seq(I0, max_age, by = step)
  B <- numeric(length(ages))
  B[1] <- B0
  for (k in 2:length(ages)) {
    B[k] <- project_B(ages[k - 1], B[k - 1], ages[k], S, alpha)
  }
  prog <- data.frame(site = S, age = ages,
                     basal_area_m2ha = B,
                     volume_m3ha = static_volume(ages, S, B, beta))
  prog <- calculate_ica_ima(prog)
  prog$extrapolated <- prog$age > max_obs_age
  prog
}

# ICA(I) = V(I) - V(I - 1)   (needs a row exactly 1 year earlier -> NA before that)
# IMA(I) = V(I) / I
# The age step must divide 1 year (e.g. 1, 0.5, 0.25) so that "I - 1" exists.
calculate_ica_ima <- function(prog) {
  step <- min(diff(prog$age))
  if (abs(1 / step - round(1 / step)) > 1e-8) {
    stop("AGE_STEP must divide 1 year (1, 0.5, 0.25, ...) so ICA = V(I) - V(I-1) exists.")
  }
  idx <- vapply(prog$age - 1, function(x) {
    m <- which(abs(prog$age - x) < 1e-8)
    if (length(m) == 0) NA_integer_ else m[1]
  }, integer(1))
  prog$ica_m3ha <- prog$volume_m3ha - prog$volume_m3ha[idx]
  prog$ima_m3ha <- prog$volume_m3ha / prog$age
  prog
}

# Technical rotation age = age where ICA crosses IMA (ICA(I) = IMA(I)).
# Equation 1 has the property that projecting I0 -> I1 -> I2 in two steps gives
# the same B2 as projecting I0 -> I2 in one step, so V(I) can be evaluated at
# ANY age directly from (I0, B0). That lets us find the crossing exactly
# (uniroot) instead of only on the simulation grid.
# Returns NA if no crossing is found before age_limit.
technical_rotation_age <- function(S, I0, B0, alpha, beta, age_limit = 40) {
  V_at <- function(I) project_V(I0, B0, I, S, alpha, beta)
  diff_fn  <- function(I) (V_at(I) - V_at(I - 1)) - V_at(I) / I      # ICA - IMA
  age_grid <- seq(I0 + 1.01, age_limit, by = 0.05)
  fv <- vapply(age_grid, diff_fn, numeric(1))
  k  <- which(fv[-length(fv)] * fv[-1] < 0)                       # sign change
  if (length(k) == 0) return(NA_real_)
  uniroot(diff_fn, c(age_grid[k[1]], age_grid[k[1] + 1]))$root
}
