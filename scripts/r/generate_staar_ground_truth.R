#!/usr/bin/env Rscript
# Regenerates the `spa_binary_full` and `ai_staar` sections of
# src/staar/testdata/ground_truth.json by running upstream STAAR:
#   - STAAR::STAAR_Binary_SPA
#   - STAAR::AI_STAAR
#
# Requires R 4.4 with Rcpp, RcppArmadillo, Matrix, GMMAT, GENESIS, STAAR,
# jsonlite. Upstream sources live at /n/netscratch/xlin/Lab/vineet/upstream.
# The STAAR stack is installed at /n/netscratch/xlin/Lab/vineet/Rlib-staar
# (separate from the user's Bioc 3.22 Rlib to avoid downgrading deps).
#
# Usage:
#   module load R/4.4.3
#   export R_LIBS_USER=/n/netscratch/xlin/Lab/vineet/Rlib-staar
#   srun -p xlin -t 00:30:00 --mem 8G \
#     Rscript scripts/r/generate_staar_ground_truth.R

suppressPackageStartupMessages({
  library(STAAR)
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
repo_root <- if (length(file_arg) > 0) {
  normalizePath(file.path(dirname(sub("^--file=", "", file_arg)), "../.."))
} else if (nzchar(Sys.getenv("FAVOR_REPO_ROOT"))) {
  Sys.getenv("FAVOR_REPO_ROOT")
} else {
  getwd()
}
json_path <- file.path(repo_root, "src", "staar", "testdata", "ground_truth.json")
stopifnot(file.exists(json_path))

mat_to_rows <- function(m) lapply(seq_len(nrow(m)), function(i) as.numeric(m[i, ]))

generate_spa_binary <- function() {
  set.seed(20261001)
  n <- 500
  p <- 10
  q <- 3

  X <- cbind(
    intercept = rep(1, n),
    age = scale(rnorm(n, mean = 55, sd = 10))[, 1],
    sex = rbinom(n, 1, 0.5)
  )
  logit_mu <- -2.6 + 0.2 * X[, "age"] + 0.05 * X[, "sex"]
  y <- rbinom(n, 1, 1 / (1 + exp(-logit_mu)))
  stopifnot(mean(y) > 0.04, mean(y) < 0.20)

  # R filters variants whose sample MAF exceeds the rare cutoff; Rust does
  # not. Oversample then keep the first p that pass so the post-filter set
  # is deterministic on both sides.
  mafs_true <- runif(5 * p, min = 0.001, max = 0.008)
  keep_cols <- integer(0)
  G_all <- matrix(0, nrow = n, ncol = length(mafs_true))
  for (j in seq_along(mafs_true)) {
    G_all[, j] <- rbinom(n, size = 2, prob = mafs_true[j])
    sample_maf <- sum(G_all[, j]) / (2 * n)
    if (sample_maf > 0 && sample_maf < 0.01) keep_cols <- c(keep_cols, j)
    if (length(keep_cols) == p) break
  }
  stopifnot(length(keep_cols) == p)
  G <- G_all[, keep_cols, drop = FALSE]

  annotation_phred <- matrix(runif(p * q, min = 0, max = 30), nrow = p, ncol = q)
  colnames(annotation_phred) <- paste0("ann", seq_len(q))

  obj_nullmodel <- fit_null_glm_Binary_SPA(y ~ age + sex,
                                           data = data.frame(y = y, age = X[, "age"], sex = X[, "sex"]),
                                           family = binomial(link = "logit"))

  res <- STAAR_Binary_SPA(G, obj_nullmodel, annotation_phred,
                          rare_maf_cutoff = 0.01, rv_num_cutoff = 2)
  stopifnot(all(res$RV_label))

  list(
    n_samples        = n,
    n_variants       = p,
    n_annotations    = q,
    case_rate        = mean(y),
    X                = mat_to_rows(X),
    y                = as.numeric(y),
    G                = mat_to_rows(G),
    annotation_phred = mat_to_rows(annotation_phred),
    mu               = as.numeric(obj_nullmodel$fitted.values),
    residuals        = as.numeric(y - obj_nullmodel$fitted.values),
    expected = list(
      num_variant   = res$num_variant,
      cMAC          = res$cMAC,
      RV_label      = as.logical(res$RV_label),
      STAAR_B       = as.numeric(res$results_STAAR_B),
      STAAR_B_1_25  = as.numeric(unlist(res$results_STAAR_B_1_25)),
      STAAR_B_1_1   = as.numeric(unlist(res$results_STAAR_B_1_1))
    )
  )
}

generate_ai_staar <- function() {
  set.seed(20261002)
  n       <- 300
  p       <- 8
  q       <- 3
  n_pops  <- 3
  B_extra <- 4
  pop_seed <- 7590

  pop   <- sort(rep(seq_len(n_pops), length.out = n)) - 1L
  pop_r <- pop + 1L

  X <- cbind(
    intercept = rep(1, n),
    age = scale(rnorm(n, mean = 50, sd = 12))[, 1],
    pc1 = scale(rnorm(n))[, 1]
  )
  y <- 2 * X[, "age"] + 0.5 * X[, "pc1"] + rnorm(n, sd = 1)

  mafs_true <- runif(5 * p, min = 0.002, max = 0.008)
  keep_cols <- integer(0)
  G_all <- matrix(0, nrow = n, ncol = length(mafs_true))
  for (j in seq_along(mafs_true)) {
    G_all[, j] <- rbinom(n, size = 2, prob = mafs_true[j])
    sample_maf <- sum(G_all[, j]) / (2 * n)
    if (sample_maf > 0 && sample_maf < 0.01) keep_cols <- c(keep_cols, j)
    if (length(keep_cols) == p) break
  }
  stopifnot(length(keep_cols) == p)
  G <- G_all[, keep_cols, drop = FALSE]

  annotation_phred <- matrix(runif(p * q, min = 0, max = 30), nrow = p, ncol = q)
  colnames(annotation_phred) <- paste0("ann", seq_len(q))

  obj_nullmodel <- fit_null_glm(y ~ age + pc1,
                                data = data.frame(y = y, age = X[, "age"], pc1 = X[, "pc1"]),
                                family = gaussian())

  # Column 1 is all-ones (the MAF-only base run); subsequent columns are
  # |N(0,1)|. Matches STAARpipeline::staar2aistaar_nullmodel and AncestryInfo.
  set.seed(pop_seed)
  make_weights <- function() {
    rnd <- abs(matrix(rnorm(n_pops * B_extra), nrow = n_pops, ncol = B_extra))
    cbind(matrix(1, nrow = n_pops, ncol = 1), rnd)
  }
  pop_weights_1_1  <- make_weights()
  pop_weights_1_25 <- make_weights()

  obj_nullmodel$pop_weights_1_1  <- pop_weights_1_1
  obj_nullmodel$pop_weights_1_25 <- pop_weights_1_25
  obj_nullmodel$pop.groups       <- pop_r

  res <- AI_STAAR(G, obj_nullmodel, annotation_phred,
                  rare_maf_cutoff = 0.01, rv_num_cutoff = 2)
  stopifnot(all(res$RV_label))

  list(
    n_samples          = n,
    n_variants         = p,
    n_annotations      = q,
    n_populations      = n_pops,
    B                  = ncol(pop_weights_1_1),
    pop_seed           = pop_seed,
    X                  = mat_to_rows(X),
    y                  = as.numeric(y),
    G                  = mat_to_rows(G),
    annotation_phred   = mat_to_rows(annotation_phred),
    pop_groups_0_based = as.integer(pop),
    pop_weights_1_1    = mat_to_rows(pop_weights_1_1),
    pop_weights_1_25   = mat_to_rows(pop_weights_1_25),
    expected = list(
      num_variant    = res$num_variant,
      cMAC           = res$cMAC,
      RV_label       = as.logical(res$RV_label),
      STAAR_O        = as.numeric(res$results_STAAR_O),
      ACAT_O         = as.numeric(res$results_ACAT_O),
      STAAR_S_1_25   = as.numeric(unlist(res$results_STAAR_S_1_25)),
      STAAR_S_1_1    = as.numeric(unlist(res$results_STAAR_S_1_1)),
      STAAR_B_1_25   = as.numeric(unlist(res$results_STAAR_B_1_25)),
      STAAR_B_1_1    = as.numeric(unlist(res$results_STAAR_B_1_1)),
      STAAR_A_1_25   = as.numeric(unlist(res$results_STAAR_A_1_25)),
      STAAR_A_1_1    = as.numeric(unlist(res$results_STAAR_A_1_1))
    )
  )
}

cat("Generating SPA binary fixture...\n")
spa_binary_full <- generate_spa_binary()

cat("Generating AI-STAAR fixture...\n")
ai_staar <- generate_ai_staar()

cat("Merging into", json_path, "\n")

# Preserve byte-identical serialization of untouched sections so the
# bit-level invariance test against staar_continuous stays stable. Only
# the two keys we regenerate get re-serialized; everything else is sliced
# out of the original file as raw text.
render_section <- function(value) {
  chunk <- jsonlite::toJSON(list(x = value), digits = 22, auto_unbox = TRUE,
                            null = "null", na = "null", pretty = 2)
  chunk <- sub('^\\s*\\{\\s*\\"x\\":\\s*', "", chunk)
  chunk <- sub('\\s*\\}\\s*$', "", chunk)
  chunk
}

raw <- readLines(json_path, warn = FALSE)
text <- paste(raw, collapse = "\n")

replace_top_level_key <- function(text, key, new_body) {
  m <- regexpr(paste0('\\"', key, '\\":'), text)
  if (m[1] < 0) stop(sprintf("key %s not found", key))
  start <- m[1] + attr(m, "match.length") - 1L
  rest <- substr(text, start + 1L, nchar(text))
  first_char <- regmatches(rest, regexpr("[^[:space:]]", rest))
  opener <- first_char
  closer <- switch(opener, "{" = "}", "[" = "]", stop("unexpected opener: ", opener))
  depth <- 0L
  in_str <- FALSE
  esc <- FALSE
  n <- nchar(rest)
  end_rel <- NA_integer_
  for (i in seq_len(n)) {
    ch <- substr(rest, i, i)
    if (esc) { esc <- FALSE; next }
    if (ch == "\\") { esc <- TRUE; next }
    if (ch == "\"") { in_str <- !in_str; next }
    if (in_str) next
    if (ch == opener) depth <- depth + 1L
    if (ch == closer) {
      depth <- depth - 1L
      if (depth == 0L) { end_rel <- i; break }
    }
  }
  stopifnot(!is.na(end_rel))
  paste0(substr(text, 1L, start),
         " ", new_body,
         substr(text, start + end_rel + 1L, nchar(text)))
}

text <- replace_top_level_key(text, "spa_binary_full", render_section(spa_binary_full))
text <- replace_top_level_key(text, "ai_staar",        render_section(ai_staar))

writeLines(text, con = json_path, sep = "")

cat("spa_binary_full: n=", spa_binary_full$n_samples,
    " case_rate=", round(spa_binary_full$case_rate, 3),
    " num_variant=", spa_binary_full$expected$num_variant, "\n", sep = "")
cat("ai_staar:        n=", ai_staar$n_samples,
    " n_pops=", ai_staar$n_populations,
    " num_variant=", ai_staar$expected$num_variant,
    " STAAR_O=", signif(ai_staar$expected$STAAR_O, 3), "\n", sep = "")
