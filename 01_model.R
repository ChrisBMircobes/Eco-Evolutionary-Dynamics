# =========================================================================
# 01_model.R
#
# Simulation engine for eco-evolutionary dynamics in a competitive
# community. Species compete under stochastic Leslie-Gower dynamics; each
# species is a pool of clones, and new clones arise by mutation in growth
# rate (r), carrying capacity (K) or competitive interactions (alpha).
#
# Competitive loads are computed from per-species aggregates, so no
# clone x clone matrix is ever built. The interaction modifier between two
# clones is exp(alpha_out_j[s_i] - alpha_in_i[s_j]), which is separable and
# always positive, so the aggregation is exact.
#
# Mutation effects for r, K and alpha share one truncated normal
# distribution, hard-capped at +/- max_effect:
#   r, K  : trait is multiplied by (1 + delta)
#   alpha : effect is added on the log scale as log1p(delta), so the
#           interaction modifier changes by exactly (1 + delta)
#
# Used by 02_run_batch.R.
# =========================================================================


# -------------------------------------------------------------------------
# Interaction matrix
# Diagonal is 1; off-diagonal entries are normal draws truncated at 0.
# -------------------------------------------------------------------------
generate_alpha_matrix <- function(nspp, mean_alpha, alpha_sd, rng_seed = NULL) {
  if (!is.null(rng_seed)) set.seed(rng_seed)
  alpha <- matrix(0, nrow = nspp, ncol = nspp)
  for (i in 1:nspp) {
    for (j in 1:nspp) {
      if (i == j) {
        alpha[i, j] <- 1.0
      } else {
        alpha[i, j] <- max(0, rnorm(1, mean = mean_alpha, sd = alpha_sd))
      }
    }
  }
  alpha
}

# Per-species summaries of the interaction matrix.
# alpha0[i, j] is the effect of species j on species i.
compute_species_traits <- function(alpha0, nspp) {
  out_strength <- vapply(1:nspp, function(i) mean(alpha0[-i, i]), numeric(1))  # effect of i on others
  in_strength  <- vapply(1:nspp, function(i) mean(alpha0[i, -i]), numeric(1))  # effect of others on i
  data.frame(
    species = 1:nspp,
    out_strength = out_strength,
    in_strength = in_strength,
    net_competitiveness = out_strength - in_strength
  )
}


# -------------------------------------------------------------------------
# Clone initialisation: one ancestral clone per species
# -------------------------------------------------------------------------
init_clones <- function(nspp, N0, K0, r0) {
  clones <- vector("list", nspp)
  for (i in 1:nspp) {
    clones[[i]] <- list(
      id = paste0("sp", i, "_ancestor"), species = i, N = as.numeric(N0[i]),
      r = r0[i], K = K0[i],
      alpha_out = rep(0, nspp), alpha_in = rep(0, nspp),
      birth_gen = 0, parent_id = "ancestor", mutation_type = "ancestor",
      delta_r = 0, delta_K = 0,
      own_delta_r = 0, own_delta_K = 0, own_alpha_effect = 0
    )
  }
  clones
}


# -------------------------------------------------------------------------
# Competitive load on each clone.
#
# The pairwise definition is
#   load_i = sum_j alpha0[s_i, s_j] * exp(out_j[s_i] - in_i[s_j]) * N_j
# Grouping clones j by species k gives
#   load_i = sum_k alpha0[s_i, k] * exp(-in_i[k]) * M[k, s_i]
#   M[k, s] = sum_{j in species k} N_j * exp(out_j[s])
# Own-species effects are never mutated (they stay 0), so intraspecific
# competition is exactly 1 per individual.
# -------------------------------------------------------------------------
compute_loads <- function(clones, alpha0) {
  n <- length(clones)
  nspp <- ncol(alpha0)
  S  <- vapply(clones, function(c) c$species, integer(1))
  Ns <- vapply(clones, function(c) c$N, numeric(1))
  AlphaOut <- matrix(unlist(lapply(clones, `[[`, "alpha_out")), nrow = n, byrow = TRUE)
  AlphaIn  <- matrix(unlist(lapply(clones, `[[`, "alpha_in")),  nrow = n, byrow = TRUE)

  tmp <- rowsum(Ns * exp(AlphaOut), S)            # (species present) x nspp
  M <- matrix(0, nspp, nspp)
  M[as.integer(rownames(tmp)), ] <- tmp           # M[k, s]

  BaseSp <- alpha0[S, , drop = FALSE]             # [i, k] = alpha0[s_i, k]
  Mt     <- t(M)[S, , drop = FALSE]               # [i, k] = M[k, s_i]
  rowSums(BaseSp * exp(-AlphaIn) * Mt)
}

# Compares the aggregated loads with a brute-force double loop on random
# clones. Returns the maximum relative error (expect ~1e-15).
check_aggregation <- function(nspp = 5, n = 60, seed = 1) {
  set.seed(seed)
  alpha0 <- generate_alpha_matrix(nspp, 0.5, 0.2)
  clones <- lapply(seq_len(n), function(i) {
    s <- sample.int(nspp, 1)
    ao <- rnorm(nspp, 0, 0.1); ai <- rnorm(nspp, 0, 0.1)
    ao[s] <- 0; ai[s] <- 0
    list(species = s, N = runif(1, 1, 100), alpha_out = ao, alpha_in = ai)
  })
  brute <- vapply(seq_len(n), function(i) {
    ci <- clones[[i]]
    sum(vapply(seq_len(n), function(j) {
      cj <- clones[[j]]
      alpha0[ci$species, cj$species] * exp(cj$alpha_out[ci$species] - ci$alpha_in[cj$species]) * cj$N
    }, numeric(1)))
  }, numeric(1))
  max(abs(brute - compute_loads(clones, alpha0)) / brute)
}


# -------------------------------------------------------------------------
# One generation of Leslie-Gower dynamics with demographic stochasticity.
# Expected N is N * (1 + r) / (1 + load / K); the realised count is Poisson
# (or Bernoulli when the expectation is below 1).
# -------------------------------------------------------------------------
update_clones_LG <- function(clones, alpha0) {
  n <- length(clones)
  if (n == 0) return(list())

  loads <- compute_loads(clones, alpha0)
  Ns <- vapply(clones, function(c) c$N, numeric(1))
  r  <- vapply(clones, function(c) c$r, numeric(1))
  K  <- vapply(clones, function(c) c$K, numeric(1))

  expected <- Ns * (1 + r) / (1 + loads / K)

  new_N <- numeric(n)
  ok    <- is.finite(expected) & expected > 0
  small <- ok & expected < 1
  big   <- ok & !small
  if (any(small)) new_N[small] <- rbinom(sum(small), 1, expected[small])
  if (any(big))   new_N[big]   <- rpois(sum(big), pmin(expected[big], 1e7))

  keep <- which(new_N > 0)
  out <- clones[keep]
  for (m in seq_along(keep)) out[[m]]$N <- new_N[keep[m]]
  out
}


# -------------------------------------------------------------------------
# Mutation effect size: normal truncated to [-max_effect, +max_effect] with
# standard deviation `sd` and mean mean_shift * sd (a negative shift makes
# deleterious mutations more likely than beneficial ones). Sampled by
# inverse CDF, so there is no spike at the cap.
# -------------------------------------------------------------------------
draw_effect <- function(n, max_effect, sd, mean_shift) {
  if (n == 0) return(numeric(0))
  m  <- mean_shift * sd
  lo <- pnorm(-max_effect, m, sd)
  hi <- pnorm( max_effect, m, sd)
  e  <- qnorm(runif(n, lo, hi), m, sd)
  pmin(max_effect, pmax(-max_effect, e))
}


# -------------------------------------------------------------------------
# Mutation step, applied to all clones once per generation.
#
# mu_vec is per species (length nspp); mut_params is
# list(max_effect, sd, mean_shift). Each mutation hits r, K or alpha with
# probability 1/3. Alpha mutations are "aggressive" (change the effect of
# the mutant on other species) or "defensive" (change the effect of others on
# the mutant); defensive_ratio is the ratio of defensive to aggressive
# mutations, so the default of 1 makes the two equally likely.
#
# The own_* fields record the effect of the mutation itself, not the
# cumulative effect inherited along the lineage (that is delta_r, delta_K).
#
# If the clone pool exceeds max_clones, the clones with the smallest
# population are dropped (ties broken at random).
# -------------------------------------------------------------------------
apply_mutations <- function(clones, generation, clone_counter, mu_vec, nspp, mut_params,
                             max_clones = 20000, K_ceiling = Inf, defensive_ratio = 1) {
  n <- length(clones)
  if (n == 0) return(list(clones = clones, clone_counter = clone_counter, peak = 0L, capped = FALSE))

  S  <- vapply(clones, function(c) c$species, integer(1))
  Ns <- vapply(clones, function(c) c$N, numeric(1))
  n_mut_vec <- rpois(n, mu_vec[S] * Ns)
  total_mut <- sum(n_mut_vec)
  if (total_mut == 0) return(list(clones = clones, clone_counter = clone_counter, peak = n, capped = FALSE))

  parent_idx <- rep(seq_len(n), times = n_mut_vec)
  channel    <- sample(c("r", "K", "alpha"), total_mut, replace = TRUE)
  p_def <- defensive_ratio / (1 + defensive_ratio)
  alpha_type <- sample(c("aggressive", "defensive"), total_mut, replace = TRUE,
                       prob = c(1 - p_def, p_def))
  eff <- draw_effect(total_mut, mut_params$max_effect, mut_params$sd, mut_params$mean_shift)

  new_clones <- vector("list", total_mut)
  for (k in seq_len(total_mut)) {
    parent <- clones[[parent_idx[k]]]
    clone_counter <- clone_counter + 1
    mutant <- parent
    mutant$id <- paste0("sp", parent$species, "_mut", clone_counter)
    mutant$N <- 1
    mutant$birth_gen <- generation
    mutant$parent_id <- parent$id
    mutant$own_delta_r <- 0; mutant$own_delta_K <- 0; mutant$own_alpha_effect <- 0

    e <- eff[k]
    ch <- channel[k]
    if (ch == "r") {
      mutant$r <- max(0.001, parent$r * (1 + e))
      mutant$delta_r <- parent$delta_r + e
      mutant$own_delta_r <- e
      mutant$mutation_type <- "Evo_r"
    } else if (ch == "K") {
      mutant$K <- min(K_ceiling, max(1, parent$K * (1 + e)))
      mutant$delta_K <- parent$delta_K + e
      mutant$own_delta_K <- e
      mutant$mutation_type <- "Evo_K"
    } else {
      other_species <- setdiff(1:nspp, parent$species)
      mutant$own_alpha_effect <- e
      if (alpha_type[k] == "aggressive") {
        mutant$alpha_out[other_species] <- parent$alpha_out[other_species] + log1p(e)
        mutant$mutation_type <- "Niche_shift_aggressive"
      } else {
        mutant$alpha_in[other_species] <- parent$alpha_in[other_species] + log1p(e)
        mutant$mutation_type <- "Niche_shift_defensive"
      }
    }
    new_clones[[k]] <- mutant
  }

  all_clones <- c(clones, new_clones)
  peak <- length(all_clones)
  capped <- FALSE
  if (peak > max_clones) {
    Ns_all <- vapply(all_clones, function(c) c$N, numeric(1))
    # N is integer valued, so adding U(0, 0.5) only breaks ties
    keep_idx <- order(Ns_all + runif(peak, 0, 0.5), decreasing = TRUE)[seq_len(max_clones)]
    all_clones <- all_clones[keep_idx]
    capped <- TRUE
  }
  list(clones = all_clones, clone_counter = clone_counter, peak = peak, capped = capped)
}


# -------------------------------------------------------------------------
# Diversity indices from a vector of species frequencies
# -------------------------------------------------------------------------
calculate_diversity <- function(frequencies) {
  alive <- frequencies[frequencies > 0]
  if (length(alive) == 0) return(c(richness = 0, shannon = 0, simpson = 0, evenness = 0))
  richness <- length(alive)
  shannon  <- -sum(alive * log(alive))
  simpson  <- sum(alive^2)
  evenness <- if (richness > 1) shannon / log(richness) else NA
  c(richness = richness, shannon = shannon, simpson = simpson, evenness = evenness)
}


# -------------------------------------------------------------------------
# Selection coefficient per species: the slope of ln(N) over the first
# window_gens generations of a run without mutations, minus the mean slope
# across species.
# -------------------------------------------------------------------------
compute_selection_coefficients <- function(community_states_no_mut, nspp, window_gens = 50) {
  df <- community_states_no_mut[community_states_no_mut$gen <= window_gens, ]
  slopes <- vapply(1:nspp, function(sp) {
    sub <- df[df$species == sp & df$N > 0, ]
    if (nrow(sub) < 2 || length(unique(sub$gen)) < 2) return(NA_real_)
    fit <- tryCatch(lm(log(N) ~ gen, data = sub), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    unname(coef(fit)["gen"])
  }, numeric(1))
  slopes - mean(slopes, na.rm = TRUE)
}


# -------------------------------------------------------------------------
# Full simulation run.
#
# Species totals are recorded every record_interval generations and every
# generation up to early_record_gens. If track_mutation_freq is TRUE, each
# mutant lineage that reaches mutation_freq_threshold of its species'
# population at any point is logged once.
#
# Returns a list with:
#   community_states   gen, species, N, freq at recorded generations
#   final_clones       clone list at the end of the run
#   n_generations_run  last generation simulated
#   mutation_log       one row per lineage that crossed the threshold
#   peak_clones        largest clone pool seen (before any truncation)
#   cap_hits           number of generations in which max_clones was binding
# -------------------------------------------------------------------------
run_simulation <- function(nspp, gens, N0, K0, r0, alpha0, with_mutations, mu_vec, mut_params,
                            record_interval = 25, max_clones = 20000, early_record_gens = 50,
                            K_ceiling = NULL, track_mutation_freq = TRUE,
                            mutation_freq_threshold = 0.10) {

  if (is.null(K_ceiling)) K_ceiling <- 10 * max(K0)

  clones <- init_clones(nspp, N0, K0, r0)
  clone_counter <- nspp
  community_states <- vector("list", 0)

  crossed_ids <- character(0)
  crossing_records <- vector("list", 0)
  peak_clones <- 0L
  cap_hits <- 0L

  t_final <- 0
  for (t in seq_len(gens)) {
    t_final <- t
    if (length(clones) == 0) break
    if (t %% 100 == 0) gc(verbose = FALSE)

    clones <- update_clones_LG(clones, alpha0)

    if (with_mutations && length(clones) > 0) {
      res <- apply_mutations(clones, t, clone_counter, mu_vec, nspp, mut_params, max_clones, K_ceiling)
      clones <- res$clones
      clone_counter <- res$clone_counter
      peak_clones <- max(peak_clones, res$peak)
      if (res$capped) cap_hits <- cap_hits + 1L
    }
    if (length(clones) == 0) break
    peak_clones <- max(peak_clones, length(clones))

    Sc  <- vapply(clones, function(c) c$species, integer(1))
    Ncl <- vapply(clones, function(c) c$N, numeric(1))
    tmp <- rowsum(Ncl, Sc)
    species_totals <- numeric(nspp)
    species_totals[as.integer(rownames(tmp))] <- tmp[, 1]
    total <- sum(species_totals)
    species_freqs <- if (total > 0) species_totals / total else rep(0, nspp)

    if ((t %% record_interval == 0 || t <= early_record_gens) && total > 0) {
      community_states[[length(community_states) + 1]] <- data.frame(
        gen = t, species = 1:nspp, N = species_totals, freq = species_freqs
      )
    }

    # Log mutant lineages the first time they exceed the within-species
    # frequency threshold
    if (with_mutations && track_mutation_freq) {
      is_mut <- vapply(clones, function(c) c$mutation_type != "ancestor", logical(1))
      sp_tot <- species_totals[Sc]
      clone_freq <- ifelse(sp_tot > 0, Ncl / sp_tot, 0)
      newly_cross <- is_mut & clone_freq >= mutation_freq_threshold
      if (any(newly_cross)) {
        ids <- vapply(clones, function(c) c$id, character(1))
        not_yet <- newly_cross & !(ids %in% crossed_ids)
        if (any(not_yet)) {
          idx <- which(not_yet)
          new_ids <- ids[idx]
          crossed_ids <- c(crossed_ids, new_ids)
          sub_clones <- clones[idx]
          crossing_records[[length(crossing_records) + 1]] <- data.frame(
            mutant_id = new_ids,
            species = vapply(sub_clones, function(c) c$species, integer(1)),
            parent_id = vapply(sub_clones, function(c) c$parent_id, character(1)),
            mutation_type = vapply(sub_clones, function(c) c$mutation_type, character(1)),
            own_delta_r = vapply(sub_clones, function(c) c$own_delta_r, numeric(1)),
            own_delta_K = vapply(sub_clones, function(c) c$own_delta_K, numeric(1)),
            own_alpha_effect = vapply(sub_clones, function(c) c$own_alpha_effect, numeric(1)),
            birth_gen = vapply(sub_clones, function(c) c$birth_gen, numeric(1)),
            gen_crossed_10pct = t,
            freq_at_crossing = clone_freq[idx],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  community_df <- if (length(community_states) > 0) do.call(rbind, community_states) else data.frame()

  # Combine the threshold-crossing records with each lineage's final state
  mutation_log <- data.frame()
  if (with_mutations && track_mutation_freq && length(crossing_records) > 0) {
    crossing_df <- do.call(rbind, crossing_records)

    if (length(clones) > 0) {
      final_ids <- vapply(clones, function(c) c$id, character(1))
      final_N   <- vapply(clones, function(c) c$N, numeric(1))
      final_species <- vapply(clones, function(c) c$species, integer(1))
      final_species_totals <- numeric(nspp)
      for (s in 1:nspp) final_species_totals[s] <- sum(final_N[final_species == s])
    } else {
      final_ids <- character(0); final_N <- numeric(0); final_species <- integer(0)
      final_species_totals <- numeric(nspp)
    }

    lookup_idx <- match(crossing_df$mutant_id, final_ids)
    alive_at_end <- !is.na(lookup_idx)
    fN <- ifelse(alive_at_end, final_N[lookup_idx], 0)
    fS <- ifelse(alive_at_end, final_species[lookup_idx], NA_integer_)
    fFreq <- ifelse(alive_at_end & final_species_totals[fS] > 0, fN / final_species_totals[fS], 0)

    outcome_df <- data.frame(
      mutant_id = crossing_df$mutant_id,
      final_gen = t_final, final_N = fN, final_freq = fFreq,
      alive_at_end = alive_at_end, stringsAsFactors = FALSE
    )
    mutation_log <- merge(crossing_df, outcome_df, by = "mutant_id", all.x = TRUE)
  }

  list(community_states = community_df, final_clones = clones, n_generations_run = t_final,
       mutation_log = mutation_log, peak_clones = peak_clones, cap_hits = cap_hits)
}


# -------------------------------------------------------------------------
# Final per-species state (N, frequency, extinct flag) from a finished run
# -------------------------------------------------------------------------
extract_final_state <- function(sim_result, nspp) {
  cd <- sim_result$community_states
  if (nrow(cd) == 0) return(data.frame(species = 1:nspp, N = 0, freq = 0, extinct = TRUE))
  final_gen <- max(cd$gen)
  final <- cd[cd$gen == final_gen, ]
  data.frame(species = final$species, N = final$N, freq = final$freq, extinct = final$N == 0)
}
