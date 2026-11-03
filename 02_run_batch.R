# =========================================================================
# 02_run_batch.R
#
# Runs one task of a batch of eco-evolutionary simulations (model in
# 01_model.R). For each parameter set, the community is simulated
# N_STOCH_REPS times with mutations and N_STOCH_REPS times without, and the
# outcomes are summarised per species.
#
# Parameter sets
#   Task `task_id` runs sets (task_id - 1) * REPLICATES_PER_TASK + 1 ...
#   task_id * REPLICATES_PER_TASK. If DESIGN_FILE is set (a table made by
#   00_make_design.R), mean_alpha, alpha_cv and k_decades are read from the
#   matching row; otherwise they are drawn at random from their ranges.
#   Per-species K, r and mu are always drawn at random within each set. With
#   a design file, the K range is centred on the geometric centre of
#   K_MIN..K_MAX and has a width of k_decades (log10 units).
#
# Initial abundance
#   N0_VAL gives every species the same starting abundance. If N0_FRAC > 0,
#   each species instead starts at N0_FRAC * r0 * K0 and N0_VAL is ignored.
#
# Outputs (written to the working directory)
#   species_data_task<id>.csv     one row per species per parameter set
#   raw_freq_traj_task<id>.rds    per-repetition species frequencies at
#                                 each checkpoint and at the end
#   mutation_log_task<id>.csv     mutant lineages that reached
#                                 MUTATION_FREQ_THRESHOLD in their species
#
# Usage: Rscript 02_run_batch.R <task_id>
#   All settings are read from environment variables (see below); the
#   task id falls back to SLURM_ARRAY_TASK_ID if not given as an argument.
# =========================================================================

source("01_model.R")

args <- commandArgs(trailingOnly = TRUE)
task_id <- if (length(args) >= 1) as.integer(args[1]) else as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "1"))

NSPP <- as.numeric(Sys.getenv("NSPP", "10"))
N0_VAL <- as.numeric(Sys.getenv("N0_VAL", "50000"))
N0_FRAC <- as.numeric(Sys.getenv("N0_FRAC", "0"))      # > 0: N0 = N0_FRAC * r0 * K0 per species (overrides N0_VAL)
GENS <- as.numeric(Sys.getenv("GENS", "1000"))
RECORD_INTERVAL <- as.numeric(Sys.getenv("RECORD_INTERVAL", "25"))
EARLY_WINDOW <- as.numeric(Sys.getenv("EARLY_WINDOW", "50"))
MAX_CLONES <- as.numeric(Sys.getenv("MAX_CLONES", "2000000"))
K_CEILING_MULT <- as.numeric(Sys.getenv("K_CEILING_MULT", "100"))
REPLICATES_PER_TASK <- as.numeric(Sys.getenv("REPLICATES_PER_TASK", "1"))
BASE_SEED <- as.numeric(Sys.getenv("BASE_SEED", "300000"))
N_STOCH_REPS <- as.numeric(Sys.getenv("N_STOCH_REPS", "100"))
CHECKPOINT_INTERVAL <- as.numeric(Sys.getenv("CHECKPOINT_INTERVAL", "50"))  # must be a multiple of RECORD_INTERVAL
MUTATION_FREQ_THRESHOLD <- as.numeric(Sys.getenv("MUTATION_FREQ_THRESHOLD", "0.10"))
NULL_TEST <- as.logical(Sys.getenv("NULL_TEST", "FALSE"))  # TRUE: the "with mutations" runs have no mutations (null baseline)

# ---- mutation effect distribution -------------------------------------
MAX_EFFECT <- as.numeric(Sys.getenv("MAX_EFFECT", "0.3"))
mut_sd_env <- Sys.getenv("MUT_SD", "")
MUT_SD <- if (nzchar(mut_sd_env)) as.numeric(mut_sd_env) else MAX_EFFECT / 3
MUT_MEAN_SHIFT <- as.numeric(Sys.getenv("MUT_MEAN_SHIFT", "0"))
MUT_PARAMS <- list(max_effect = MAX_EFFECT, sd = MUT_SD, mean_shift = MUT_MEAN_SHIFT)

# ---- parameter ranges ---------------------------------------------------
K_MIN <- as.numeric(Sys.getenv("K_MIN", "1e5")); K_MAX <- as.numeric(Sys.getenv("K_MAX", "1e6"))
R_MIN <- as.numeric(Sys.getenv("R_MIN", "0.05")); R_MAX <- as.numeric(Sys.getenv("R_MAX", "0.2"))
MU_MIN <- as.numeric(Sys.getenv("MU_MIN", "1e-5")); MU_MAX <- as.numeric(Sys.getenv("MU_MAX", "1e-3"))
ALPHA_MEAN_MIN <- as.numeric(Sys.getenv("ALPHA_MEAN_MIN", "0.05")); ALPHA_MEAN_MAX <- as.numeric(Sys.getenv("ALPHA_MEAN_MAX", "1.0"))
ALPHA_CV_MIN <- as.numeric(Sys.getenv("ALPHA_CV_MIN", "0.1")); ALPHA_CV_MAX <- as.numeric(Sys.getenv("ALPHA_CV_MAX", "0.5"))

# ---- optional design table ---------------------------------------------
DESIGN_FILE <- Sys.getenv("DESIGN_FILE", "")
design <- NULL
if (nzchar(DESIGN_FILE)) {
  design <- read.csv(DESIGN_FILE)
  missing_cols <- setdiff(c("mean_alpha", "alpha_cv", "k_decades"), names(design))
  if (length(missing_cols) > 0) stop("DESIGN_FILE is missing columns: ", paste(missing_cols, collapse = ", "))
}

if (RECORD_INTERVAL > 0 && CHECKPOINT_INTERVAL %% RECORD_INTERVAL != 0) {
  stop("CHECKPOINT_INTERVAL must be a multiple of RECORD_INTERVAL, or the checkpoint generations won't be recorded")
}
CHECKPOINTS <- seq(CHECKPOINT_INTERVAL, GENS, by = CHECKPOINT_INTERVAL)
QUANTILE_PROBS <- c(0.10, 0.25, 0.50, 0.75, 0.90)

# Column-wise quantiles of a (reps x species) matrix, returned as a data
# frame with one row per species
build_quantile_df <- function(mat, prefix) {
  qmat <- apply(mat, 2, quantile, probs = QUANTILE_PROBS, na.rm = TRUE)
  qmat <- t(qmat)
  colnames(qmat) <- paste0(prefix, "_p", QUANTILE_PROBS * 100)
  as.data.frame(qmat)
}

# (N_STOCH_REPS x NSPP) frequency matrix -> long-format per-repetition rows
freq_matrix_to_block <- function(mat, gen_label, scenario, nspp) {
  df <- as.data.frame(mat)
  colnames(df) <- paste0("freq_sp", 1:nspp)
  df$rep_i <- seq_len(nrow(mat))
  df$gen <- gen_label
  df$scenario <- scenario
  df
}

# One checkpoint slice of a (reps x species x checkpoints) array, kept as a
# matrix even when there is only one repetition (R would otherwise drop
# the dimensions)
slice_ck <- function(arr, ck) matrix(arr[, , ck], nrow = dim(arr)[1])

cat("=== Task", task_id, "running", REPLICATES_PER_TASK, "replicates x", N_STOCH_REPS,
    "reps, checkpoints at", paste(CHECKPOINTS, collapse = ","), "===\n")
cat("Mutation effects: truncated normal, cap +/-", MAX_EFFECT, ", sd", MUT_SD,
    ", mean shift", MUT_MEAN_SHIFT, "SD; MAX_CLONES", MAX_CLONES, "\n")
if (is.null(design)) {
  cat("Alpha: mean", ALPHA_MEAN_MIN, "-", ALPHA_MEAN_MAX, ", CV", ALPHA_CV_MIN, "-", ALPHA_CV_MAX,
      " (random draws, no design file)\n")
} else {
  cat("Design file", DESIGN_FILE, "with", nrow(design), "rows: mean_alpha, alpha_cv, k_decades from the table\n")
}
cat("N0:", if (N0_FRAC > 0) paste0(N0_FRAC, " x r0*K0 per species") else paste0(N0_VAL, " per species"),
    "; NULL_TEST =", NULL_TEST, "\n")

all_rows <- vector("list", REPLICATES_PER_TASK)
all_traj_blocks <- vector("list", REPLICATES_PER_TASK)
all_mut_logs <- vector("list", REPLICATES_PER_TASK)

for (local_i in seq_len(REPLICATES_PER_TASK)) {
  replicate_id <- (task_id - 1) * REPLICATES_PER_TASK + local_i
  seed <- BASE_SEED + replicate_id
  set.seed(seed)

  # K range for this set: from the design (width in decades, centred on the
  # K_MIN..K_MAX centre) or the full K_MIN..K_MAX range
  if (!is.null(design)) {
    if (replicate_id > nrow(design)) stop("replicate_id ", replicate_id, " exceeds the ", nrow(design), " rows of ", DESIGN_FILE)
    drow <- design[replicate_id, ]
    k_decades <- drow$k_decades
    k_centre  <- (log10(K_MIN) + log10(K_MAX)) / 2
    K_lo <- 10^(k_centre - k_decades / 2); K_hi <- 10^(k_centre + k_decades / 2)
  } else {
    k_decades <- log10(K_MAX / K_MIN)
    K_lo <- K_MIN; K_hi <- K_MAX
  }

  K0 <- exp(runif(NSPP, log(K_lo), log(K_hi)))
  r0 <- runif(NSPP, R_MIN, R_MAX)
  mu_vec <- exp(runif(NSPP, log(MU_MIN), log(MU_MAX)))
  if (!is.null(design)) {
    mean_alpha <- drow$mean_alpha
    alpha_cv   <- drow$alpha_cv
  } else {
    mean_alpha <- runif(1, ALPHA_MEAN_MIN, ALPHA_MEAN_MAX)
    alpha_cv   <- runif(1, ALPHA_CV_MIN, ALPHA_CV_MAX)
  }
  alpha_sd   <- alpha_cv * mean_alpha
  alpha0 <- generate_alpha_matrix(NSPP, mean_alpha, alpha_sd)
  realized_mean_alpha <- mean(alpha0[row(alpha0) != col(alpha0)])
  species_traits <- compute_species_traits(alpha0, NSPP)
  K_ceiling <- K_CEILING_MULT * max(K0)

  N0 <- if (N0_FRAC > 0) N0_FRAC * r0 * K0 else rep(N0_VAL, NSPP)

  outcome_counts <- matrix(0L, nrow = NSPP, ncol = 4,
                            dimnames = list(NULL, c("extinct_both", "rescued", "evolution_caused_loss", "survives_both")))
  n_increased_freq <- integer(NSPP)
  sel_coef_mat <- matrix(NA_real_, nrow = N_STOCH_REPS, ncol = NSPP)
  delta_shannon_vec <- numeric(N_STOCH_REPS)
  richness_no_vec <- numeric(N_STOCH_REPS)
  richness_with_vec <- numeric(N_STOCH_REPS)
  peak_clones_vec <- numeric(N_STOCH_REPS)
  cap_hits_vec <- numeric(N_STOCH_REPS)

  freq_final_no   <- matrix(NA_real_, N_STOCH_REPS, NSPP)
  freq_final_with <- matrix(NA_real_, N_STOCH_REPS, NSPP)
  freq_ckpt_no   <- array(NA_real_, dim = c(N_STOCH_REPS, NSPP, length(CHECKPOINTS)))
  freq_ckpt_with <- array(NA_real_, dim = c(N_STOCH_REPS, NSPP, length(CHECKPOINTS)))
  richness_final_no_int  <- integer(N_STOCH_REPS)
  richness_final_with_int <- integer(N_STOCH_REPS)

  rep_mut_logs <- vector("list", N_STOCH_REPS)

  for (rep_i in seq_len(N_STOCH_REPS)) {
    dyn_seed <- BASE_SEED + replicate_id * 10000L + rep_i

    # "with mutations" run (or, if NULL_TEST, a no-mutation run with a different seed)
    set.seed(dyn_seed + if (NULL_TEST) 1L else 0L)
    sim_with <- run_simulation(nspp = NSPP, gens = GENS, N0 = N0, K0 = K0, r0 = r0, alpha0 = alpha0,
                                with_mutations = !NULL_TEST, mu_vec = mu_vec, mut_params = MUT_PARAMS,
                                record_interval = RECORD_INTERVAL, max_clones = MAX_CLONES,
                                early_record_gens = EARLY_WINDOW, K_ceiling = K_ceiling,
                                track_mutation_freq = TRUE,
                                mutation_freq_threshold = MUTATION_FREQ_THRESHOLD)
    set.seed(dyn_seed)
    sim_no <- run_simulation(nspp = NSPP, gens = GENS, N0 = N0, K0 = K0, r0 = r0, alpha0 = alpha0,
                              with_mutations = FALSE, mu_vec = mu_vec, mut_params = MUT_PARAMS,
                              record_interval = RECORD_INTERVAL, max_clones = MAX_CLONES,
                              early_record_gens = EARLY_WINDOW, K_ceiling = K_ceiling,
                              track_mutation_freq = FALSE)

    final_with <- extract_final_state(sim_with, NSPP)
    final_no   <- extract_final_state(sim_no, NSPP)
    sel_coef_mat[rep_i, ] <- compute_selection_coefficients(sim_no$community_states, NSPP, window_gens = EARLY_WINDOW)
    peak_clones_vec[rep_i] <- sim_with$peak_clones
    cap_hits_vec[rep_i]    <- sim_with$cap_hits

    div_with <- calculate_diversity(final_with$freq)
    div_no   <- calculate_diversity(final_no$freq)
    delta_shannon_vec[rep_i] <- unname(div_with["shannon"] - div_no["shannon"])
    richness_no_vec[rep_i]   <- unname(div_no["richness"])
    richness_with_vec[rep_i] <- unname(div_with["richness"])

    freq_final_no[rep_i, ] <- final_no$freq
    freq_final_with[rep_i, ] <- final_with$freq
    richness_final_no_int[rep_i]  <- sum(!final_no$extinct)
    richness_final_with_int[rep_i] <- sum(!final_with$extinct)

    cs_no <- sim_no$community_states; cs_with <- sim_with$community_states
    for (ck in seq_along(CHECKPOINTS)) {
      g <- CHECKPOINTS[ck]
      row_no <- cs_no[cs_no$gen == g, ]
      if (nrow(row_no) == NSPP) freq_ckpt_no[rep_i, , ck] <- row_no$freq[order(row_no$species)]
      row_with <- cs_with[cs_with$gen == g, ]
      if (nrow(row_with) == NSPP) freq_ckpt_with[rep_i, , ck] <- row_with$freq[order(row_with$species)]
    }

    # Per-species outcome category, comparing runs without and with mutations
    for (s in 1:NSPP) {
      en <- final_no$extinct[s]; ew <- final_with$extinct[s]
      cat_s <- if (en && ew) "extinct_both"
               else if (en && !ew) "rescued"
               else if (!en && ew) "evolution_caused_loss"
               else "survives_both"
      outcome_counts[s, cat_s] <- outcome_counts[s, cat_s] + 1L
      if (final_with$freq[s] > final_no$freq[s]) n_increased_freq[s] <- n_increased_freq[s] + 1L
    }

    if (nrow(sim_with$mutation_log) > 0) {
      ml <- sim_with$mutation_log
      ml$task_id <- task_id; ml$replicate_id <- replicate_id; ml$rep_i <- rep_i
      rep_mut_logs[[rep_i]] <- ml
    }
  }

  # ---- quantile columns (species level, NSPP rows) -----------------------
  qdf_list <- list(
    build_quantile_df(freq_final_no,   "freq_final_no"),
    build_quantile_df(freq_final_with, "freq_final_with")
  )
  for (ck in seq_along(CHECKPOINTS)) {
    g <- CHECKPOINTS[ck]
    qdf_list[[length(qdf_list) + 1]] <- build_quantile_df(slice_ck(freq_ckpt_no, ck),   paste0("freq_gen", g, "_no"))
    qdf_list[[length(qdf_list) + 1]] <- build_quantile_df(slice_ck(freq_ckpt_with, ck), paste0("freq_gen", g, "_with"))
  }
  quantile_cols_df <- do.call(cbind, qdf_list)

  # ---- richness histogram (one value per set, repeated on every species row)
  richness_hist_no  <- sapply(0:NSPP, function(k) sum(richness_final_no_int == k))
  richness_hist_with <- sapply(0:NSPP, function(k) sum(richness_final_with_int == k))
  names(richness_hist_no) <- paste0("n_survivors_", 0:NSPP, "_no")
  names(richness_hist_with) <- paste0("n_survivors_", 0:NSPP, "_with")
  richness_hist_df <- as.data.frame(as.list(c(richness_hist_no, richness_hist_with)))

  rows <- data.frame(
    task_id = task_id, replicate_id = replicate_id, species = 1:NSPP,
    K0 = K0, r0 = r0, mu = mu_vec, N0 = N0,
    mutation_supply = K0 * mu_vec,
    net_competitiveness = species_traits$net_competitiveness,
    out_strength = species_traits$out_strength, in_strength = species_traits$in_strength,
    mean_alpha = mean_alpha, alpha_cv = alpha_cv, alpha_sd = alpha_sd,
    realized_mean_alpha = realized_mean_alpha,
    k_decades = k_decades, K_min_eff = K_lo, K_max_eff = K_hi,
    design_used = !is.null(design),
    max_effect = MAX_EFFECT, mut_sd = MUT_SD, mut_mean_shift = MUT_MEAN_SHIFT,
    null_test = NULL_TEST,
    n_reps = N_STOCH_REPS,
    n_extinct_both = outcome_counts[, "extinct_both"],
    n_rescued = outcome_counts[, "rescued"],
    n_evolution_caused_loss = outcome_counts[, "evolution_caused_loss"],
    n_survives_both = outcome_counts[, "survives_both"],
    n_increased_freq = n_increased_freq,
    selection_coef_mean = colMeans(sel_coef_mat, na.rm = TRUE),
    selection_coef_sd   = apply(sel_coef_mat, 2, sd, na.rm = TRUE),
    mean_delta_shannon = mean(delta_shannon_vec), sd_delta_shannon = sd(delta_shannon_vec),
    mean_richness_no = mean(richness_no_vec), mean_richness_with = mean(richness_with_vec),
    peak_clones_max = max(peak_clones_vec),
    cap_hits_total = sum(cap_hits_vec)
  )
  rows <- cbind(rows, quantile_cols_df)
  rows <- cbind(rows, richness_hist_df[rep(1, NSPP), ])
  all_rows[[local_i]] <- rows

  # ---- raw per-repetition trajectory blocks
  traj_blocks <- vector("list", 2 * (length(CHECKPOINTS) + 1))
  bi <- 0
  for (ck in seq_along(CHECKPOINTS)) {
    bi <- bi + 1; traj_blocks[[bi]] <- freq_matrix_to_block(slice_ck(freq_ckpt_no, ck),   CHECKPOINTS[ck], "no",   NSPP)
    bi <- bi + 1; traj_blocks[[bi]] <- freq_matrix_to_block(slice_ck(freq_ckpt_with, ck), CHECKPOINTS[ck], "with", NSPP)
  }
  bi <- bi + 1; traj_blocks[[bi]] <- freq_matrix_to_block(freq_final_no,   "final", "no",   NSPP)
  bi <- bi + 1; traj_blocks[[bi]] <- freq_matrix_to_block(freq_final_with, "final", "with", NSPP)
  traj_df <- do.call(rbind, traj_blocks)
  traj_df$task_id <- task_id
  traj_df$replicate_id <- replicate_id
  all_traj_blocks[[local_i]] <- traj_df

  # ---- mutation log for this parameter set
  rep_mut_logs <- rep_mut_logs[!vapply(rep_mut_logs, is.null, logical(1))]
  if (length(rep_mut_logs) > 0) all_mut_logs[[local_i]] <- do.call(rbind, rep_mut_logs)

  cat("  replicate", local_i, "/", REPLICATES_PER_TASK, "done; mean alpha", round(mean_alpha, 3),
      "(realized", round(realized_mean_alpha, 3), "); K range", format(K_lo, digits = 3), "-",
      format(K_hi, digits = 3), "; peak clones", max(peak_clones_vec),
      "; cap hits", sum(cap_hits_vec), "\n")
}

out <- do.call(rbind, all_rows)
write.csv(out, paste0("species_data_task", task_id, ".csv"), row.names = FALSE)
cat("Wrote species_data_task", task_id, ".csv (", nrow(out), "rows,", ncol(out), "columns )\n")

traj_out <- do.call(rbind, all_traj_blocks)
saveRDS(traj_out, paste0("raw_freq_traj_task", task_id, ".rds"), compress = "xz")
cat("Wrote raw_freq_traj_task", task_id, ".rds (", nrow(traj_out), "rows )\n")

mut_logs_nonnull <- all_mut_logs[!vapply(all_mut_logs, is.null, logical(1))]
mut_out <- if (length(mut_logs_nonnull) > 0) do.call(rbind, mut_logs_nonnull) else data.frame()
write.csv(mut_out, paste0("mutation_log_task", task_id, ".csv"), row.names = FALSE)
cat("Wrote mutation_log_task", task_id, ".csv (", nrow(mut_out), "rows )\n")
