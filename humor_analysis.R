# =====================================================================
# Re-run of all QUALITATIVE statistics after recoding
# Inputs (same folder):
#   CODING  : your recoded sheet (one row per translation; two coder columns with comma-separated labels)
#   SCORES  : "JOKES - Sheet5.csv"   (rater scores)
# Packages: lme4, blme (psych optional). No `irr` needed - alpha is implemented below.
# Edit ONLY the CONFIG block and the three dictionaries.
# =====================================================================
suppressPackageStartupMessages({library(lme4); library(blme)})
set.seed(42)

# ------------------------------- CONFIG -------------------------------
CODING <- "Data - Data.csv"          # your recoded file
SCORES <- "JOKES - Sheet5.csv"
COL_CODER1 <- 3; COL_CODER2 <- 4; COL_MODEL <- 5     # column POSITIONS in the coding sheet
# Direction/joke are inferred from row order, as in the original sheet:
# rows 1-36 = Urdu->English (6 models x 6 jokes), rows 37-72 = English->Urdu; jokes cycle 1..6 within a model.
# If your recoded sheet has explicit Direction / Joke columns, set these to their positions instead:
COL_DIR <- NA; COL_JOKE <- NA
closed <- c("gpt-5.2-2025-12-11","gemini2.5-flash","claude-sonnet-4-5")

# ---------------------------- DICTIONARIES ----------------------------
# 1) raw label (lower-case, bracket tag removed)  ->  canonical theme name
alias <- c("inaccuracy"="Content Modification", "language"="Code-Mixing",
           "lack of cultural context"="Cultural Context", "literality lose humour"="Literality Losing Humor",
           "literality loses humor"="Literality Losing Humor", "literality losing humour"="Literality Losing Humor",
           "flow"="Flow and Naturalness", "flow and naturalness"="Flow and Naturalness",
           "respect"="Respect and Formality", "respect and formality"="Respect and Formality",
           "incoherency"="Incoherency","content modification"="Content Modification","lexical choice"="Lexical Choice",
           "cultural context"="Cultural Context","structural addition"="Structural Addition",
           "grammatical errors"="Grammatical Errors","code-mixing"="Code-Mixing")
# 2) meta-theme of each canonical theme (EDIT to match your recoded Tables 5/6)
meta <- c("Incoherency"="Mechanical","Content Modification"="Mechanical","Grammatical Errors"="Mechanical",
          "Code-Mixing"="Mechanical","Respect and Formality"="Mechanical",
          "Lexical Choice"="Conceptual","Literality Losing Humor"="Conceptual","Flow and Naturalness"="Conceptual",
          "Cultural Context"="Conceptual","Structural Addition"="Conceptual")
# 3) polarity when a label has NO [Positive]/[Negative] tag (Structural Addition is always positive in this corpus)
default_pol <- c("Structural Addition"="Positive")        # confirmed: every Structural Addition code is positive. Everything else defaults to "Negative"
# Alternative taxonomy for the sensitivity analysis (Content Modification treated as a meaning error)
meta_alt <- meta; meta_alt["Content Modification"] <- "Conceptual"

# ------------------------------ PARSING -------------------------------
cod <- read.csv(CODING, stringsAsFactors = FALSE, check.names = FALSE)
sc  <- read.csv(SCORES, stringsAsFactors = FALSE)
n <- nrow(cod); stopifnot(n == 72)
if (is.na(COL_DIR)) { cod$Direction <- rep(c("Urdu->English","English->Urdu"), each = 36); cod$Joke.ID <- rep(1:6, length.out = n)
} else { cod$Direction <- cod[[COL_DIR]]; cod$Joke.ID <- cod[[COL_JOKE]] }
cod$Model <- cod[[COL_MODEL]]
cod$class <- ifelse(cod$Model %in% closed, "closed", "open")

UNTAGGED <- character()
parse_cell <- function(x) {                     # "Flow [Positive], Lexical Choice" -> data.frame(theme, pol)
  if (is.na(x) || trimws(x) == "") return(data.frame(theme = character(), pol = character()))
  p <- trimws(strsplit(x, ",")[[1]]); p <- p[p != ""]
  pol <- ifelse(grepl("posit", p, ignore.case = TRUE), "Positive", ifelse(grepl("negat", p, ignore.case = TRUE), "Negative", NA))
  raw <- tolower(trimws(gsub("\\[.*$", "", p)))
  th <- unname(alias[raw])
  bad <- is.na(th); if (any(bad)) warning("UNMAPPED label(s): ", paste(unique(p[bad]), collapse = " | "))
  pol <- ifelse(is.na(pol), ifelse(th %in% names(default_pol), default_pol[th], "Negative"), pol)
  # flag bipolar themes that were left untagged
  untag <- is.na(ifelse(grepl("posit|negat", p, ignore.case = TRUE), 1, NA)) & th %in% c("Structural Addition","Cultural Context","Lexical Choice","Flow and Naturalness")
  if (any(untag[!bad])) UNTAGGED <<- c(UNTAGGED, th[untag])
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
cat("\nLabels with NO [Positive]/[Negative] tag (default polarity applied) - tag these in your sheet if they should be positive:\n"); print(table(UNTAGGED))

# ---------------------- 1. INTER-CODER RELIABILITY ---------------------
# Krippendorff's alpha for two coders, complete data, arbitrary distance d(a,b).
kalpha <- function(A, B, dist) {                 # A, B: lists of values (one per unit)
  N <- length(A); vals <- c(A, B); M <- length(vals)
  Do <- mean(mapply(dist, A, B))
  tot <- 0; for (i in 1:M) for (j in 1:M) if (i != j) tot <- tot + dist(vals[[i]], vals[[j]])
  1 - Do / (tot / (M * (M - 1)))
}
d_nominal <- function(a, b) as.numeric(!identical(a, b))
kalpha_bin <- function(a, b) {                    # fast closed form for two coders, binary 0/1 data
  vals <- c(a, b); M <- length(vals); n1 <- sum(vals)
  Do <- mean(a != b); De <- 2 * n1 * (M - n1) / (M * (M - 1)); if (De == 0) NA else 1 - Do / De }
d_masi <- function(a, b) { if (length(a) == 0 && length(b) == 0) return(0)
  i <- length(intersect(a, b)); u <- length(union(a, b)); j <- i / u
  m <- if (setequal(a, b)) 1 else if (i == length(a) || i == length(b)) 2/3 else if (i > 0) 1/3 else 0
  1 - j * m }
wide <- function(v) { x <- reshape(long[, c("unit","coder","theme", v)], idvar = c("unit","theme"), timevar = "coder", direction = "wide"); x }
w <- wide("present"); names(w)[3:4] <- c("c1","c2")
# (a) per theme: binary, nominal distance + Cohen's kappa + % agreement
rel <- do.call(rbind, lapply(themes, function(t) { s <- subset(w, theme == t)
po <- mean(s$c1 == s$c2); pe <- mean(s$c1) * mean(s$c2) + (1 - mean(s$c1)) * (1 - mean(s$c2))
data.frame(theme = t, n_c1 = sum(s$c1), n_c2 = sum(s$c2), pct_agree = po,
           kappa = if (pe == 1) NA else (po - pe) / (1 - pe),
           alpha_nominal = kalpha_bin(s$c1, s$c2)) }))
cat("\n== 1a. Per-theme reliability ==\n"); print(rel, digits = 3, row.names = FALSE)
# (b) all themes pooled: every unit x theme cell is one binary observation (what a single overall alpha usually means)
cat("\n== 1b. Pooled binary alpha (all unit x theme cells):", round(kalpha_bin(w$c1, w$c2), 3), "\n")
# (c) set-based: each translation's whole LABEL SET is one unit, MASI distance (respects multi-label coding)
sets <- lapply(1:2, function(cd) lapply(seq_len(n), function(i) themes[long$present[long$unit == i & long$coder == cd] == 1]))
cat("== 1c. Set-based alpha (MASI distance, unit = translation):", round(kalpha(sets[[1]], sets[[2]], d_masi), 3), "\n")
# (d) bootstrap CI for the pooled alpha, resampling TRANSLATIONS
idx <- split(seq_len(nrow(w)), w$unit)
bs <- replicate(2000, { r <- unlist(idx[sample(n, replace = TRUE)]); kalpha_bin(w$c1[r], w$c2[r]) })
cat("   pooled alpha 95% bootstrap CI:", round(quantile(bs, c(.025, .975)), 3), "\n")
write.csv(rel, "out_reliability_by_theme.csv", row.names = FALSE)

# ---------- 2. FINAL (CONSENSUS = UNION) CODING PER TRANSLATION ----------
u <- aggregate(cbind(present, pos, neg) ~ unit + theme, long, max)
tw <- reshape(u[, c("unit","theme","present")], idvar = "unit", timevar = "theme", direction = "wide"); names(tw) <- sub("present\\.", "", names(tw))
pw <- reshape(u[, c("unit","theme","pos")], idvar = "unit", timevar = "theme", direction = "wide"); names(pw) <- sub("pos\\.", "pos_", names(pw))
fin <- merge(merge(cbind(cod[, c("Direction","Joke.ID","Model","class")], unit = 1:n), tw, by = "unit"), pw, by = "unit")
fin$Direction <- factor(fin$Direction, levels = c("English->Urdu","Urdu->English"))
fin$any_pos <- as.integer(rowSums(fin[paste0("pos_", themes)]) > 0)           # any positive-polarity code
fin$MECH <- as.integer(rowSums(fin[names(meta)[meta == "Mechanical"]]) > 0)
fin$CONC <- as.integer(rowSums(fin[names(meta)[meta == "Conceptual"]]) > 0)
fin$MECH_alt <- as.integer(rowSums(fin[names(meta_alt)[meta_alt == "Mechanical"]]) > 0)
fin$CONC_alt <- as.integer(rowSums(fin[names(meta_alt)[meta_alt == "Conceptual"]]) > 0)
for (k in c("Accuracy","Fluency","Humour","Equiv")) sc[[paste0(k, "_mean")]] <- (sc[[k]] + sc[[paste0(k, "1")]]) / 2
fin <- merge(fin, sc[, c("Model","Direction","Joke.ID", paste0(c("Accuracy","Fluency","Humour","Equiv"), "_mean"))],
             by = c("Model","Direction","Joke.ID")); stopifnot(nrow(fin) == 72)
fin$joke <- paste(fin$Direction, fin$Joke.ID, sep = "_J"); fin$class <- factor(fin$class, levels = c("open","closed"))

# ----------------------- 3. THEME COUNT TABLE -----------------------
cnt <- do.call(rbind, lapply(c(themes, "any_pos", "MECH", "CONC"), function(t) {
  col <- if (t %in% themes) t else t
  data.frame(theme = t, meta = ifelse(t %in% themes, meta[t], "composite"),
             closed_ENUR = sum(fin[[col]][fin$class == "closed" & fin$Direction == "English->Urdu"]),
             open_ENUR   = sum(fin[[col]][fin$class == "open"   & fin$Direction == "English->Urdu"]),
             closed_URENG = sum(fin[[col]][fin$class == "closed" & fin$Direction == "Urdu->English"]),
             open_URENG   = sum(fin[[col]][fin$class == "open"   & fin$Direction == "Urdu->English"])) }))
cat("\n== 3. Counts out of 18 translations per cell (union of coders) ==\n"); print(cnt, row.names = FALSE)
cat("\nPositive-polarity counts (translations with a [Positive] code), per theme and cell:\n")
pos_tab <- do.call(rbind, lapply(themes, function(t) data.frame(theme = t,
                                                                closed_ENUR = sum(fin[[paste0("pos_", t)]][fin$class == "closed" & fin$Direction == "English->Urdu"]),
                                                                open_ENUR   = sum(fin[[paste0("pos_", t)]][fin$class == "open"   & fin$Direction == "English->Urdu"]),
                                                                closed_URENG = sum(fin[[paste0("pos_", t)]][fin$class == "closed" & fin$Direction == "Urdu->English"]),
                                                                open_URENG   = sum(fin[[paste0("pos_", t)]][fin$class == "open"   & fin$Direction == "Urdu->English"]))))
print(pos_tab, row.names = FALSE); write.csv(pos_tab, "out_positive_counts.csv", row.names = FALSE)
write.csv(cnt, "out_theme_counts.csv", row.names = FALSE)

# ------------------ 4. CO-OCCURRENCE OF MECHANICAL / CONCEPTUAL ------------------
cat("\n== 4. Mechanical x Conceptual co-occurrence (your taxonomy) ==\n")
print(table(class = fin$class, Direction = fin$Direction, Mechanical = fin$MECH, Conceptual = fin$CONC))
cat("\nSensitivity: Content Modification moved to Conceptual\n")
print(table(class = fin$class, Direction = fin$Direction, Mechanical = fin$MECH_alt, Conceptual = fin$CONC_alt))

# --------------------- 5. THEME x SCORE LINKS (descriptive) ---------------------
sc_cols <- paste0(c("Accuracy","Fluency","Humour","Equiv"), "_mean")
link <- do.call(rbind, lapply(c(themes, "any_pos", "MECH", "CONC"), function(t) do.call(rbind, lapply(c("closed","open"), function(cl) {
  s <- fin[fin$class == cl, ]; a <- s[s[[t]] == 1, ]; b <- s[s[[t]] == 0, ]
  if (nrow(a) < 2 || nrow(b) < 2) return(NULL)
  data.frame(theme = t, class = cl, n_with = nrow(a), n_without = nrow(b),
             setNames(as.list(round(colMeans(a[sc_cols]) - colMeans(b[sc_cols]), 2)), paste0("diff_", sub("_mean","",sc_cols)))) }))))
cat("\n== 5. Mean score DIFFERENCE (with theme - without), within class ==\n"); print(link, row.names = FALSE)
write.csv(link, "out_theme_score_links.csv", row.names = FALSE)

# -------- 6. PENALISED MIXED LOGISTIC: class effect and class x direction --------
targets <- c(themes, "any_pos", "MECH", "CONC")
glm_one <- function(th) {
  x <- fin; x$y <- x[[th]]
  if (sum(x$y) < 4 || sum(1 - x$y) < 4) return(NULL)                    # too rare / too common to model
  g <- function(rhs) bglmer(as.formula(paste("y ~", rhs, "+ (1|joke)")), data = x, family = binomial, fixef.prior = normal(sd = 3))
  m0 <- g("Direction"); m1 <- g("Direction + class"); m2 <- g("Direction * class")
  l <- function(a, b) pchisq(2 * (as.numeric(logLik(b)) - as.numeric(logLik(a))), 1, lower.tail = FALSE)
  data.frame(theme = th, closed_n = sum(x$y[x$class == "closed"]), open_n = sum(x$y[x$class == "open"]),
             class_p = l(m0, m1), interaction_p = l(m1, m2),
             singular = any(isSingular(m0), isSingular(m1), isSingular(m2))) }      # TRUE = joke variance estimated at ~0 (boundary fit); model is then ~ a penalised plain logistic
tt <- do.call(rbind, lapply(targets, function(th) tryCatch(glm_one(th), error = function(e) { message(th, " failed: ", conditionMessage(e)); NULL })))
tt$class_p_BH <- p.adjust(tt$class_p, "BH"); tt$interaction_p_BH <- p.adjust(tt$interaction_p, "BH")
cat("\n== 6. Theme ~ class x direction (penalised mixed logistic; LRT approximate) ==\n"); print(tt, digits = 3, row.names = FALSE)
cat("Themes not modelled (fewer than 4 occurrences or non-occurrences):", paste(setdiff(targets, tt$theme), collapse = ", "), "\n")
write.csv(tt, "out_theme_glmm.csv", row.names = FALSE)