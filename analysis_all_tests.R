# =============================================================================
# Evaluating Large Language Models for Translating Humor in English and Urdu
# Complete statistical analysis - every test, table and figure statistic in the paper
#
# Inputs (put both files in the working directory):
#   JOKES - Sheet5.csv   one row per translation (72 rows). Columns used:
#                        Direction, Joke ID, Model, and the two raters' scores:
#                        Humour/Humour1, Fluency/Fluency1, Accuracy/Accuracy1, Equiv/Equiv1
#   Data - Data.csv      qualitative coding sheet (72 rows; coder 1 in column 3,
#                        coder 2 in column 4, model in column 5)
#
# Outputs: printed to the console and written as CSV files into ./results/
#
# Packages: base R only for sections 1-7; lme4 + blme for section 8.
#   install.packages(c("lme4", "blme"))
#
# Sections (the paper table/figure each one feeds):
#   1. Load data
#   2. Descriptive scores per model, direction and metric ............ Table 5
#   3. Rater reliability: ICC, weighted kappa, Spearman, agreement,
#      rater leniency (Wilcoxon) ........................................ Results (reliability)
#   4. Friedman permutation tests + Kendall's W ......................... Table 6
#   5. Closed-source minus open-weight gap (exact Wilcoxon, bootstrap),
#      pooled and repeated for each rater ............................... Table 7
#   6. Direction tests (exact Mann-Whitney) and gap-by-direction ........ Table 8
#   7. Correlations between metrics within model class
#   8. Qualitative coding: Krippendorff's alpha, counts, co-occurrence,
#      theme-score links, mixed-effects logistic regression ............. Tables 4, 9; Figures 1-2
# =============================================================================

# ------------------------------- configuration --------------------------------
SCORES <- "JOKES - Sheet5.csv"
CODING <- "Data - Data.csv"
OUTDIR <- "results"

CLOSED <- c("gpt-5.2-2025-12-11", "gemini2.5-flash", "claude-sonnet-4-5")
OPEN   <- c("llama3.1:latest", "qwen3:8b", "deepseek-r1:latest")

B_MAIN  <- 10000   # bootstrap resamples of jokes
B_ICC   <- 5000    # bootstrap resamples of translations (reliability CIs)
B_PERM  <- 20000   # permutations for the Friedman tests
B_ALPHA <- 2000    # bootstrap resamples of translations (Krippendorff's alpha)
SEED    <- 42

METRICS <- c("humor", "fluency", "accuracy", "equivalence")
# metric -> (rater A column, rater B column) in the scores sheet
METRIC_COLS <- list(humor       = c("Humour",   "Humour1"),
                    fluency     = c("Fluency",  "Fluency1"),
                    accuracy    = c("Accuracy", "Accuracy1"),
                    equivalence = c("Equiv",    "Equiv1"))

dir.create(OUTDIR, showWarnings = FALSE)
save_csv <- function(x, name) write.csv(x, file.path(OUTDIR, name), row.names = FALSE)

# ------------------------------- helpers --------------------------------------
boot_ci <- function(x, stat = mean, B = B_MAIN) {
  n <- length(x)
  quantile(replicate(B, stat(x[sample.int(n, n, TRUE)])), c(.025, .975), names = FALSE)
}

# exact two-sided Wilcoxon signed-rank test: zero differences dropped, tied ranks averaged
exact_wsr <- function(diffs) {
  dd <- diffs[diffs != 0]; n <- length(dd)
  if (n == 0) return(c(Wplus = NA, p = NA))
  r <- rank(abs(dd)); Wobs <- sum(r[dd > 0]); half <- sum(r) / 2
  signs <- as.matrix(expand.grid(rep(list(c(0, 1)), n)))
  W <- as.vector(signs %*% r)
  c(Wplus = Wobs, p = mean(abs(W - half) >= abs(Wobs - half) - 1e-9))
}

# exact two-sided rank-sum (Mann-Whitney) test over all splits, midranks for ties
exact_ranksum <- function(x, y) {
  z <- c(x, y); r <- rank(z); n <- length(x)
  E <- n * (length(z) + 1) / 2; Robs <- sum(r[seq_len(n)])
  R <- apply(combn(length(z), n), 2, function(i) sum(r[i]))
  mean(abs(R - E) >= abs(Robs - E) - 1e-9)
}

# quadratic-weighted Cohen's kappa on the full 1-5 scale
qwk <- function(x, y, k = 5) {
  O <- table(factor(x, 1:k), factor(y, 1:k)); O <- O / sum(O)
  E <- outer(rowSums(O), colSums(O))
  W <- outer(1:k, 1:k, function(i, j) ((i - j) / (k - 1))^2)
  1 - sum(W * O) / sum(W * E)
}

# ICC from an n (translations) x k (raters) matrix, two-way ANOVA mean squares:
#   ICC(A,k): two-way random effects, absolute agreement, average of k raters
#   ICC(C,k): two-way mixed effects, consistency, average of k raters
icc_k <- function(M) {
  n <- nrow(M); k <- ncol(M); gm <- mean(M)
  MSR <- k * sum((rowMeans(M) - gm)^2) / (n - 1)
  MSC <- n * sum((colMeans(M) - gm)^2) / (k - 1)
  SSE <- sum((M - outer(rowMeans(M), rep(1, k)) - outer(rep(1, n), colMeans(M)) + gm)^2)
  MSE <- SSE / ((n - 1) * (k - 1))
  c(ICC_Ak = (MSR - MSE) / (MSR + (MSC - MSE) / n), ICC_Ck = (MSR - MSE) / MSR)
}

# =============================== 1. LOAD DATA ==================================
set.seed(SEED)
sc <- read.csv(SCORES, stringsAsFactors = FALSE)
stopifnot(nrow(sc) == 72, all(sc$Model %in% c(CLOSED, OPEN)))

raw <- sc
raw$direction <- ifelse(raw$Direction == "English->Urdu", "EN-UR", "UR-EN")
raw$joke      <- paste0(ifelse(raw$direction == "EN-UR", "EN", "UR"), raw$Joke.ID)
raw$class     <- ifelse(raw$Model %in% OPEN, "open", "closed")
raw$model     <- raw$Model

# long format: one row per translation x rater x metric
d <- do.call(rbind, lapply(METRICS, function(m) {
  cols <- METRIC_COLS[[m]]
  do.call(rbind, lapply(1:2, function(r)
    data.frame(joke = raw$joke, direction = raw$direction, model = raw$model,
               class = raw$class, rater = c("A", "B")[r], metric = m,
               score = raw[[cols[r]]], stringsAsFactors = FALSE)))
}))

# one row per translation x metric, both raters side by side
wide <- merge(subset(d, rater == "A"), subset(d, rater == "B"),
              by = c("joke", "direction", "model", "class", "metric"),
              suffixes = c(".A", ".B"))
stopifnot(nrow(wide) == 72 * length(METRICS))
wide$score.AB <- (wide$score.A + wide$score.B) / 2   # two-rater mean = the unit used for tests

# ======================= 2. DESCRIPTIVES (Table 5) =============================
# unit = joke (n = 6 per model x direction x metric); score = mean of the two raters.
# 95% t-interval truncated to the 1-5 scale; bootstrap interval shown alongside.
set.seed(SEED)
desc <- do.call(rbind, lapply(split(wide, list(wide$model, wide$direction, wide$metric), drop = TRUE),
  function(w) {
    x <- w$score.AB; n <- length(x); bs <- boot_ci(x)
    tt <- mean(x) + c(-1, 1) * qt(.975, n - 1) * sd(x) / sqrt(n)
    data.frame(model = w$model[1], class = w$class[1], direction = w$direction[1],
               metric = w$metric[1], n = n, mean = mean(x), sd = sd(x), median = median(x),
               t_lo = max(1, tt[1]), t_hi = min(5, tt[2]),
               boot_lo = bs[1], boot_hi = bs[2], row.names = NULL)
  }))
cat("\n== 2. Table 5: descriptive scores ==\n"); print(desc, digits = 3, row.names = FALSE)
save_csv(desc, "table5_descriptives.csv")

# ================== 3. RATER RELIABILITY AND LENIENCY ==========================
set.seed(SEED)
agree <- do.call(rbind, lapply(METRICS, function(m) {
  w <- subset(wide, metric == m); n <- nrow(w)
  M <- cbind(w$score.A, w$score.B)
  idx <- replicate(B_ICC, sample.int(n, n, TRUE), simplify = FALSE)
  bi <- do.call(rbind, lapply(idx, function(i) icc_k(M[i, , drop = FALSE])))
  bq <- vapply(idx, function(i) qwk(w$score.A[i], w$score.B[i]), numeric(1))
  bs <- vapply(idx, function(i) suppressWarnings(cor(w$score.A[i], w$score.B[i], method = "spearman")), numeric(1))
  ic <- icc_k(M)
  lw <- suppressWarnings(wilcox.test(w$score.A, w$score.B, paired = TRUE, exact = FALSE))
  q  <- function(v, p) quantile(v, p, names = FALSE, na.rm = TRUE)
  data.frame(metric = m, n = n,
             ICC_Ak = ic["ICC_Ak"], ICC_Ak_lo = q(bi[, "ICC_Ak"], .025), ICC_Ak_hi = q(bi[, "ICC_Ak"], .975),
             ICC_Ck = ic["ICC_Ck"], ICC_Ck_lo = q(bi[, "ICC_Ck"], .025), ICC_Ck_hi = q(bi[, "ICC_Ck"], .975),
             qwk = qwk(w$score.A, w$score.B), qwk_lo = q(bq, .025), qwk_hi = q(bq, .975),
             rho = cor(w$score.A, w$score.B, method = "spearman"), rho_lo = q(bs, .025), rho_hi = q(bs, .975),
             exact_agree = mean(w$score.A == w$score.B),
             within1 = mean(abs(w$score.A - w$score.B) <= 1),
             mean_A = mean(w$score.A), mean_B = mean(w$score.B),
             leniency_wilcoxon_p = lw$p.value, row.names = NULL)
}))
agree$leniency_p_BH <- p.adjust(agree$leniency_wilcoxon_p, "BH")
cat("\n== 3. Rater reliability (ICC(A,k), ICC(C,k), weighted kappa, Spearman, agreement, leniency) ==\n")
print(agree, digits = 3, row.names = FALSE)
save_csv(agree, "rater_reliability.csv")

# ==================== 4. FRIEDMAN PERMUTATION TESTS (Table 6) ===================
# joke = block (n = 6), six models = treatments (k = 6), per direction x metric.
# p-values: 20,000 permutations of the model labels within each joke (not the chi-square approximation).
# Kendall's W = chi-square / (n (k - 1)).
friedman_perm <- function(mat, B = B_PERM) {
  n <- nrow(mat); k <- ncol(mat)
  R <- t(apply(mat, 1, rank))                       # within-joke ranks, ties averaged
  denom <- sum(R^2) - n * k * (k + 1)^2 / 4         # tie-corrected denominator (same for every permutation)
  if (denom == 0) return(c(chisq = NA, W = NA, p = NA))
  stat <- function(R) (k - 1) * sum((colSums(R) - n * (k + 1) / 2)^2) / denom
  obs  <- stat(R)
  perm <- replicate(B, stat(t(apply(R, 1, sample))))  # shuffle model labels within each joke
  c(chisq = obs, W = obs / (n * (k - 1)), p = (1 + sum(perm >= obs - 1e-9)) / (B + 1))
}
set.seed(SEED)
fried <- do.call(rbind, lapply(c("EN-UR", "UR-EN"), function(dr)
  do.call(rbind, lapply(METRICS, function(m) {
    w <- subset(wide, metric == m & direction == dr)
    mat <- tapply(w$score.AB, list(w$joke, w$model), mean)
    stopifnot(nrow(mat) == 6, ncol(mat) == 6)
    r <- friedman_perm(mat)
    data.frame(direction = dr, metric = m, n_blocks = nrow(mat), k_models = ncol(mat),
               chisq = r["chisq"], W = r["W"], p_perm = r["p"], row.names = NULL)
  }))))
fried$p_BH <- p.adjust(fried$p_perm, "BH")          # 8 tests
cat("\n== 4. Table 6: Friedman permutation tests ==\n"); print(fried, digits = 3, row.names = FALSE)
save_csv(fried, "table6_friedman.csv")

# ============== 5. CLOSED - OPEN GAP, PAIRED BY JOKE (Table 7) ==================
# joke-level gap = mean(closed models) - mean(open models); exact Wilcoxon signed-rank on the gaps.
joke_gap <- function(dat, scorecol) {
  m <- tapply(dat[[scorecol]], list(dat$joke, dat$class), mean)
  m[, "closed"] - m[, "open"]
}
gap_row <- function(g, label, scope, metric) {
  e <- exact_wsr(g); ci <- boot_ci(g)
  data.frame(rater = label, scope = scope, metric = metric, n_jokes = length(g),
             gap = mean(g), lo = ci[1], hi = ci[2],
             closed_higher = paste0(sum(g > 0), "/", length(g)),
             Wplus = e["Wplus"], p_exact = e["p"], row.names = NULL)
}
set.seed(SEED)
gaps <- do.call(rbind, lapply(c("score.AB", "score.A", "score.B"), function(s) {
  lab <- c(score.AB = "mean of raters", score.A = "A", score.B = "B")[s]
  do.call(rbind, lapply(METRICS, function(m) {
    rbind(
      gap_row(joke_gap(subset(wide, metric == m & direction == "EN-UR"), s), lab, "EN-UR", m),
      gap_row(joke_gap(subset(wide, metric == m & direction == "UR-EN"), s), lab, "UR-EN", m),
      gap_row(joke_gap(subset(wide, metric == m), s),                        lab, "pooled", m))
  }))
}))
gaps$p_BH <- ave(gaps$p_exact, gaps$rater, FUN = function(p) p.adjust(p, "BH"))   # 12 tests per rater
cat("\n== 5. Table 7: closed-minus-open gap (rows 'mean of raters'); per-rater repeat for A and B ==\n")
print(gaps, digits = 3, row.names = FALSE)
save_csv(gaps, "table7_closed_open_gap.csv")

# ================== 6. DIRECTION TESTS (Table 8) ================================
# Different jokes in each direction, so the test is unpaired: exact Mann-Whitney over all
# C(12,6) = 924 splits of the 12 joke-level means (average of the 3 models in the class).
joke_cls_mean <- function(m, cls, dr) {
  w <- subset(wide, metric == m & class == cls & direction == dr)
  tapply(w$score.AB, w$joke, mean)
}
dir_tests <- do.call(rbind, lapply(c("closed", "open"), function(cls)
  do.call(rbind, lapply(METRICS, function(m) {
    x <- joke_cls_mean(m, cls, "EN-UR"); y <- joke_cls_mean(m, cls, "UR-EN")
    data.frame(class = cls, metric = m, n_ENUR = length(x), n_URENG = length(y),
               mean_ENUR = mean(x), mean_URENG = mean(y),
               delta_URENG_minus_ENUR = mean(y) - mean(x),
               p_exact = exact_ranksum(x, y), row.names = NULL)
  }))))
dir_tests$p_BH <- p.adjust(dir_tests$p_exact, "BH")       # 8 tests
cat("\n== 6. Table 8: direction effect within model class ==\n"); print(dir_tests, digits = 3, row.names = FALSE)
save_csv(dir_tests, "table8_direction_by_class.csv")

# does the closed-open gap differ between directions?
int_tests <- do.call(rbind, lapply(METRICS, function(m) {
  gx <- joke_gap(subset(wide, metric == m & direction == "EN-UR"), "score.AB")
  gy <- joke_gap(subset(wide, metric == m & direction == "UR-EN"), "score.AB")
  data.frame(metric = m, gap_ENUR = mean(gx), gap_URENG = mean(gy),
             diff = mean(gx) - mean(gy), p_exact = exact_ranksum(gx, gy), row.names = NULL)
}))
int_tests$p_BH <- p.adjust(int_tests$p_exact, "BH")       # 4 tests
cat("\n== 6b. Gap by direction ==\n"); print(int_tests, digits = 3, row.names = FALSE)
save_csv(int_tests, "gap_by_direction.csv")

# ================= 7. CORRELATIONS BETWEEN METRICS (within class) ===============
# Cross-rater Spearman (A's score on one metric vs B's score on the other, and the reverse),
# averaged, to avoid a single rater's halo effect. Cluster bootstrap over jokes.
metric_cor <- function(x_metric, y_metric, cls) {
  g <- function(m) { w <- subset(wide, metric == m & class == cls); w[order(w$joke, w$model), ] }
  X <- g(x_metric); Y <- g(y_metric)
  stopifnot(all(X$joke == Y$joke), all(X$model == Y$model))
  jokes <- unique(X$joke)
  stat <- function(js) {
    idx <- unlist(lapply(js, function(j) which(X$joke == j)))
    mean(c(cor(X$score.A[idx], Y$score.B[idx], method = "spearman"),
           cor(X$score.B[idx], Y$score.A[idx], method = "spearman")))
  }
  est <- stat(jokes)
  bs  <- replicate(B_MAIN, suppressWarnings(stat(sample(jokes, length(jokes), TRUE))))
  c(rho = est, lo = quantile(bs, .025, names = FALSE, na.rm = TRUE),
    hi = quantile(bs, .975, names = FALSE, na.rm = TRUE))
}
set.seed(SEED)
mcor <- do.call(rbind, lapply(c("closed", "open"), function(cls)
  do.call(rbind, lapply(list(c("accuracy", "humor"), c("equivalence", "humor"), c("fluency", "humor")),
    function(pr) {
      r <- metric_cor(pr[1], pr[2], cls)
      data.frame(class = cls, metric_x = pr[1], metric_y = pr[2],
                 rho = r[1], lo = r[2], hi = r[3], row.names = NULL)
    }))))
cat("\n== 7. Spearman correlations between metrics ==\n"); print(mcor, digits = 3, row.names = FALSE)
save_csv(mcor, "metric_correlations.csv")

# ======================= 8. QUALITATIVE CODING ================================
suppressPackageStartupMessages({library(lme4); library(blme)})
set.seed(SEED)

COL_CODER1 <- 3; COL_CODER2 <- 4; COL_MODEL <- 5     # columns in the coding sheet
COL_DIR <- NA; COL_JOKE <- NA                        # NA = rows are ordered 36 UR->EN then 36 EN->UR, jokes 1-6 repeating

# map the free-text labels in the coding sheet to the ten themes
alias <- c("inaccuracy"="Content Modification", "language"="Code-Mixing",
           "lack of cultural context"="Cultural Context", "literality lose humour"="Literality Losing Humor",
           "literality loses humor"="Literality Losing Humor", "literality losing humour"="Literality Losing Humor",
           "flow"="Flow and Naturalness", "flow and naturalness"="Flow and Naturalness",
           "respect"="Respect and Formality", "respect and formality"="Respect and Formality",
           "incoherency"="Incoherency", "content modification"="Content Modification", "lexical choice"="Lexical Choice",
           "cultural context"="Cultural Context", "structural addition"="Structural Addition",
           "grammatical errors"="Grammatical Errors", "code-mixing"="Code-Mixing")
meta <- c("Incoherency"="Mechanical", "Content Modification"="Mechanical", "Grammatical Errors"="Mechanical",
          "Code-Mixing"="Mechanical", "Respect and Formality"="Mechanical",
          "Lexical Choice"="Conceptual", "Literality Losing Humor"="Conceptual", "Flow and Naturalness"="Conceptual",
          "Cultural Context"="Conceptual", "Structural Addition"="Conceptual")
default_pol <- c("Structural Addition"="Positive")   # polarity when a label carries no [Positive]/[Negative] tag

cod <- read.csv(CODING, stringsAsFactors = FALSE, check.names = FALSE)
n <- nrow(cod); stopifnot(n == 72)
if (is.na(COL_DIR)) {
  cod$Direction <- rep(c("Urdu->English", "English->Urdu"), each = 36); cod$Joke.ID <- rep(1:6, length.out = n)
} else { cod$Direction <- cod[[COL_DIR]]; cod$Joke.ID <- cod[[COL_JOKE]] }
cod$Model <- cod[[COL_MODEL]]
cod$class <- ifelse(cod$Model %in% CLOSED, "closed", "open")

UNTAGGED <- character()
parse_cell <- function(x) {
  if (is.na(x) || trimws(x) == "") return(data.frame(theme = character(), pol = character()))
  p <- trimws(strsplit(x, ",")[[1]]); p <- p[p != ""]
  pol <- ifelse(grepl("posit", p, ignore.case = TRUE), "Positive", ifelse(grepl("negat", p, ignore.case = TRUE), "Negative", NA))
  raw_lab <- tolower(trimws(gsub("\\[.*$", "", p)))
  th <- unname(alias[raw_lab])
  bad <- is.na(th); if (any(bad)) warning("UNMAPPED label(s): ", paste(unique(p[bad]), collapse = " | "))
  pol <- ifelse(is.na(pol), ifelse(th %in% names(default_pol), default_pol[th], "Negative"), pol)
  untag <- !grepl("posit|negat", p, ignore.case = TRUE) & th %in% c("Structural Addition", "Cultural Context", "Lexical Choice", "Flow and Naturalness")
  if (any(untag[!bad])) UNTAGGED <<- c(UNTAGGED, th[untag & !bad])
  unique(data.frame(theme = th[!bad], pol = pol[!bad], stringsAsFactors = FALSE))
}
themes <- names(meta)
long <- do.call(rbind, lapply(seq_len(n), function(i) do.call(rbind, lapply(1:2, function(cd) {
  pc <- parse_cell(cod[i, if (cd == 1) COL_CODER1 else COL_CODER2])
  data.frame(unit = i, coder = cd, theme = themes,
             present = as.integer(themes %in% pc$theme),
             pos = as.integer(themes %in% pc$theme[pc$pol == "Positive"]),
             neg = as.integer(themes %in% pc$theme[pc$pol == "Negative"]),
             stringsAsFactors = FALSE) }))))
long$meta <- meta[long$theme]
cat("\nLabels with NO [Positive]/[Negative] tag (default polarity applied):\n"); print(table(UNTAGGED))

# ---- 8.1 inter-coder reliability (Krippendorff's alpha) ----
# per-theme binary alpha (Table 9), pooled binary alpha, and set-based alpha with MASI distance
kalpha_bin <- function(a, b) {
  vals <- c(a, b); M <- length(vals); n1 <- sum(vals)
  Do <- mean(a != b); De <- 2 * n1 * (M - n1) / (M * (M - 1)); if (De == 0) NA else 1 - Do / De
}
d_masi <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(0)
  i <- length(intersect(a, b)); u <- length(union(a, b)); j <- i / u
  m <- if (setequal(a, b)) 1 else if (i == length(a) || i == length(b)) 2/3 else if (i > 0) 1/3 else 0
  1 - j * m
}
w_codes <- {
  x <- reshape(long[, c("unit", "coder", "theme", "present")], idvar = c("unit", "theme"),
               timevar = "coder", direction = "wide")
  names(x)[3:4] <- c("c1", "c2"); x
}
rel <- do.call(rbind, lapply(themes, function(t) {
  s <- subset(w_codes, theme == t)
  po <- mean(s$c1 == s$c2); pe <- mean(s$c1) * mean(s$c2) + (1 - mean(s$c1)) * (1 - mean(s$c2))
  data.frame(theme = t, meta = meta[t], times_coded = sum(s$c1) + sum(s$c2),
             disagreements = sum(s$c1 != s$c2), pct_agree = po,
             kappa = if (pe == 1) NA else (po - pe) / (1 - pe),
             alpha = kalpha_bin(s$c1, s$c2), row.names = NULL)
}))
cat("\n== 8.1a Table 9: per-theme reliability ==\n"); print(rel, digits = 3, row.names = FALSE)
save_csv(rel, "table9_reliability_by_theme.csv")
cat("== 8.1b Pooled binary alpha (all unit x theme cells):", round(kalpha_bin(w_codes$c1, w_codes$c2), 3), "\n")

# set-based alpha: each translation's whole label set is one unit, MASI distance (respects multi-label coding)
sets <- lapply(1:2, function(cd) lapply(seq_len(n), function(i) themes[long$present[long$unit == i & long$coder == cd] == 1]))
vals <- c(sets[[1]], sets[[2]]); M <- length(vals)
D <- matrix(0, M, M)
for (i in seq_len(M)) for (j in seq_len(M)) D[i, j] <- d_masi(vals[[i]], vals[[j]])
alpha_masi <- function(s) {                          # s = indices of the translations used (with repeats when bootstrapping)
  idx <- c(s, s + n); Mi <- length(idx)
  Do <- mean(D[cbind(s, s + n)]); De <- sum(D[idx, idx]) / (Mi * (Mi - 1))
  1 - Do / De
}
a_masi <- alpha_masi(seq_len(n))
bs_masi <- replicate(B_ALPHA, alpha_masi(sample.int(n, n, TRUE)))
ci_masi <- quantile(bs_masi, c(.025, .975), names = FALSE)
cat("== 8.1c Overall set-based alpha (MASI distance, unit = translation):", round(a_masi, 3),
    " 95% bootstrap CI [", round(ci_masi[1], 3), ",", round(ci_masi[2], 3), "]\n")
save_csv(data.frame(alpha_MASI = a_masi, lo = ci_masi[1], hi = ci_masi[2]), "alpha_masi_overall.csv")

# ---- 8.2 final (consensus = union of both coders) coding per translation ----
u <- aggregate(cbind(present, pos, neg) ~ unit + theme, long, max)
tw <- reshape(u[, c("unit", "theme", "present")], idvar = "unit", timevar = "theme", direction = "wide")
names(tw) <- sub("present\\.", "", names(tw))
pw <- reshape(u[, c("unit", "theme", "pos")], idvar = "unit", timevar = "theme", direction = "wide")
names(pw) <- sub("pos\\.", "pos_", names(pw))
fin <- merge(merge(cbind(cod[, c("Direction", "Joke.ID", "Model", "class")], unit = 1:n), tw, by = "unit"), pw, by = "unit")
fin$Direction <- factor(fin$Direction, levels = c("English->Urdu", "Urdu->English"))
fin$any_pos <- as.integer(rowSums(fin[paste0("pos_", themes)]) > 0)                     # any positive-polarity code
fin$MECH <- as.integer(rowSums(fin[names(meta)[meta == "Mechanical"]]) > 0)             # any mechanical theme
fin$CONC <- as.integer(rowSums(fin[names(meta)[meta == "Conceptual"]]) > 0)             # any conceptual theme
sc_means <- sc[, c("Model", "Direction", "Joke.ID")]
for (m in METRICS) {
  cols <- METRIC_COLS[[m]]; sc_means[[paste0(m, "_mean")]] <- (sc[[cols[1]]] + sc[[cols[2]]]) / 2
}
fin <- merge(fin, sc_means, by = c("Model", "Direction", "Joke.ID")); stopifnot(nrow(fin) == 72)
fin$joke  <- paste(fin$Direction, fin$Joke.ID, sep = "_J")
fin$class <- factor(fin$class, levels = c("open", "closed"))

# ---- 8.3 theme counts per cell (Figure 1 / Table 4) ----
tcols <- c(themes, "any_pos", "MECH", "CONC")
cnt <- do.call(rbind, lapply(tcols, function(t) {
  data.frame(theme = t, meta = ifelse(t %in% themes, meta[t], "composite"),
             closed_ENUR  = sum(fin[[t]][fin$class == "closed" & fin$Direction == "English->Urdu"]),
             open_ENUR    = sum(fin[[t]][fin$class == "open"   & fin$Direction == "English->Urdu"]),
             closed_URENG = sum(fin[[t]][fin$class == "closed" & fin$Direction == "Urdu->English"]),
             open_URENG   = sum(fin[[t]][fin$class == "open"   & fin$Direction == "Urdu->English"]),
             closed_total = sum(fin[[t]][fin$class == "closed"]), open_total = sum(fin[[t]][fin$class == "open"]),
             row.names = NULL)
}))
cat("\n== 8.3 Translations coded with each theme (out of 18 per cell, 36 per class) ==\n"); print(cnt, row.names = FALSE)
save_csv(cnt, "theme_counts.csv")

pos_tab <- do.call(rbind, lapply(themes, function(t) {
  col <- paste0("pos_", t)
  data.frame(theme = t,
             closed_ENUR  = sum(fin[[col]][fin$class == "closed" & fin$Direction == "English->Urdu"]),
             open_ENUR    = sum(fin[[col]][fin$class == "open"   & fin$Direction == "English->Urdu"]),
             closed_URENG = sum(fin[[col]][fin$class == "closed" & fin$Direction == "Urdu->English"]),
             open_URENG   = sum(fin[[col]][fin$class == "open"   & fin$Direction == "Urdu->English"]),
             row.names = NULL)
}))
cat("\nPositive-polarity counts:\n"); print(pos_tab, row.names = FALSE)
save_csv(pos_tab, "positive_counts.csv")

# ---- 8.4 co-occurrence of mechanical and conceptual problems ----
cat("\n== 8.4 Mechanical x Conceptual co-occurrence ==\n")
print(table(class = fin$class, Direction = fin$Direction, Mechanical = fin$MECH, Conceptual = fin$CONC))

# ---- 8.5 theme x score links (descriptive; Figure 2) ----
sc_cols <- paste0(METRICS, "_mean")
link <- do.call(rbind, lapply(tcols, function(t) do.call(rbind, lapply(c("closed", "open"), function(cl) {
  s <- fin[fin$class == cl, ]; a <- s[s[[t]] == 1, ]; b <- s[s[[t]] == 0, ]
  if (nrow(a) < 2 || nrow(b) < 2) return(NULL)
  data.frame(theme = t, class = cl, n_with = nrow(a), n_without = nrow(b),
             setNames(as.list(round(colMeans(a[sc_cols]) - colMeans(b[sc_cols]), 2)),
                      paste0("diff_", sub("_mean", "", sc_cols))), row.names = NULL)
}))))
cat("\n== 8.5 Mean score difference (with theme - without), within class ==\n"); print(link, row.names = FALSE)
save_csv(link, "theme_score_links.csv")

# ---- 8.6 penalised mixed-effects logistic regression (Figure 1 p-values) ----
# theme ~ direction + class + (1 | joke); likelihood-ratio tests for class and for class x direction.
# bglmer with a normal(sd = 3) prior on fixed effects keeps estimates finite when a theme
# never occurs in one class; BH-adjusted across themes.
glm_one <- function(th) {
  x <- fin; x$y <- x[[th]]
  if (sum(x$y) < 4 || sum(1 - x$y) < 4) return(NULL)              # too rare / too common to model
  g <- function(rhs) bglmer(as.formula(paste("y ~", rhs, "+ (1|joke)")), data = x, family = binomial,
                            fixef.prior = normal(sd = 3))
  m0 <- g("Direction"); m1 <- g("Direction + class"); m2 <- g("Direction * class")
  l <- function(a, b) pchisq(2 * (as.numeric(logLik(b)) - as.numeric(logLik(a))), 1, lower.tail = FALSE)
  data.frame(theme = th, closed_n = sum(x$y[x$class == "closed"]), open_n = sum(x$y[x$class == "open"]),
             class_p = l(m0, m1), interaction_p = l(m1, m2),
             singular = any(isSingular(m0), isSingular(m1), isSingular(m2)))   # TRUE = joke variance ~0 (boundary fit)
}
tt <- do.call(rbind, lapply(tcols, function(th)
  tryCatch(glm_one(th), error = function(e) { message(th, " failed: ", conditionMessage(e)); NULL })))
tt$class_p_BH <- p.adjust(tt$class_p, "BH"); tt$interaction_p_BH <- p.adjust(tt$interaction_p, "BH")
cat("\n== 8.6 Theme ~ class x direction (penalised mixed logistic; approximate LRT) ==\n"); print(tt, digits = 3, row.names = FALSE)
cat("Themes not modelled (fewer than 4 occurrences or non-occurrences):", paste(setdiff(tcols, tt$theme), collapse = ", "), "\n")
save_csv(tt, "theme_glmm.csv")

cat("\nDone. Results written to ./", OUTDIR, "/\n", sep = "")
print(sessionInfo())
