#!/usr/bin/env Rscript
# ============================================================
# Pacific hake: find "most consistent" 3-consecutive-age blocks
# using abundance-at-age = Nsamp * composition-at-age
# and variability around a shifting (rolling) mean.
#
# Outputs (written to ./output_ageblock_sensitivity/):
#   1) block_metrics_window3.csv
#   2) block_metrics_window5.csv
#   3) best_blocks_summary.csv
#   4) SS_ready_age_bins.txt
#   5) PNG plots per fleet + window
#
# Input:
#   pacific_hake_cond_age_at_len_from_SS.csv
# ============================================================

options(stringsAsFactors = FALSE)

# ---------- user settings ----------
infile <- "pacific_hake_cond_age_at_len_from_SS.csv"
outdir <- "output_ageblock_sensitivity"
block_size <- 3
windows <- c(3, 5)   # sensitivity on shifting-mean window
# ----------------------------------

if (!file.exists(infile)) {
  stop(paste0("Input file not found: ", infile,
              "\nPlace the CSV in the working directory or change 'infile'."))
}

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dat <- read.csv(infile, check.names = FALSE)

# Basic checks
req_cols <- c("year", "fleet", "Nsamp", "total_fish")
missing_cols <- setdiff(req_cols, names(dat))
if (length(missing_cols) > 0) stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))

age_cols <- grep("^a[0-9]+$", names(dat), value = TRUE)
if (length(age_cols) < 3) stop("Did not find at least 3 age composition columns named like a1, a2, ...")

ages <- as.integer(sub("^a", "", age_cols))
ord <- order(ages)
age_cols <- age_cols[ord]
ages <- ages[ord]

fleets <- sort(unique(dat$fleet))
if (length(fleets) != 2) {
  warning(paste0("Expected 2 fleets; found ", length(fleets), ": ", paste(fleets, collapse = ", ")))
}

# ---- helper: centered rolling mean with partial windows at ends ----
roll_mean_partial <- function(x, k) {
  x <- as.numeric(x)
  n <- length(x)
  half <- floor(k / 2)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1, i - half)
    hi <- min(n, i + half)
    out[i] <- mean(x[lo:hi], na.rm = TRUE)
  }
  out
}

# ---- helper: metrics for a single fleet and window ----
calc_block_metrics <- function(dat_fleet, ages, age_cols, block_size = 3, window = 3) {
  years <- sort(unique(dat_fleet$year))

  # Build year x age abundance table
  # abundance = Nsamp * (pct / total_fish)
  # (total_fish is typically 100 if pct is in percent)
  year_age_abund <- matrix(NA_real_, nrow = length(years), ncol = length(ages),
                           dimnames = list(years, ages))

  for (j in seq_along(ages)) {
    col <- age_cols[j]
    pct <- dat_fleet[[col]]
    prop <- pct / dat_fleet$total_fish
    abund <- dat_fleet$Nsamp * prop

    # One row per year in this dataset; if multiple, sum them
    tmp <- tapply(abund, dat_fleet$year, sum, na.rm = TRUE)
    year_age_abund[names(tmp), as.character(ages[j])] <- tmp
  }

  # Candidate blocks
  starts <- ages[ages <= max(ages) - (block_size - 1)]
  out <- data.frame()

  for (s in starts) {
    block_ages <- s:(s + block_size - 1)
    y <- rowSums(year_age_abund[, as.character(block_ages), drop = FALSE], na.rm = TRUE)
    m <- roll_mean_partial(y, window)
    r <- y - m

    # Robust-ish dispersion metrics on residuals (variance around shifting mean)
    MAD <- median(abs(r - median(r, na.rm = TRUE)), na.rm = TRUE)
    IQR <- as.numeric(quantile(r, 0.75, na.rm = TRUE) - quantile(r, 0.25, na.rm = TRUE))
    SD  <- sd(r, na.rm = TRUE)
    VAR <- var(r, na.rm = TRUE)

    out <- rbind(out, data.frame(
      fleet = unique(dat_fleet$fleet),
      window = window,
      start_age = s,
      end_age = s + block_size - 1,
      ages = paste(block_ages, collapse = "-"),
      MAD = MAD,
      IQR = IQR,
      SD = SD,
      VAR = VAR,
      n_years = length(years)
    ))
  }

  out
}

# ---- helper: consensus "most consistent" block (rank-sum across metrics) ----
pick_consensus <- function(metrics_df) {
  m <- metrics_df
  # ranks: smaller = better
  m$rank_MAD <- rank(m$MAD, ties.method = "min")
  m$rank_IQR <- rank(m$IQR, ties.method = "min")
  m$rank_SD  <- rank(m$SD,  ties.method = "min")
  m$rank_sum <- m$rank_MAD + m$rank_IQR + m$rank_SD
  m <- m[order(m$rank_sum, m$SD), ]
  m[1, ]
}

# ---- compute metrics ----
all_metrics <- list()

for (fl in fleets) {
  dat_f <- dat[dat$fleet == fl, ]
  for (w in windows) {
    all_metrics[[paste0("fleet", fl, "_w", w)]] <- calc_block_metrics(dat_f, ages, age_cols,
                                                                    block_size = block_size,
                                                                    window = w)
  }
}

metrics_df <- do.call(rbind, all_metrics)

# Write per-window metrics
for (w in windows) {
  fn <- file.path(outdir, paste0("block_metrics_window", w, ".csv"))
  write.csv(metrics_df[metrics_df$window == w, ], fn, row.names = FALSE)
}

# ---- summarize "best" blocks ----
best_rows <- data.frame()

for (fl in fleets) {
  for (w in windows) {
    m <- metrics_df[metrics_df$fleet == fl & metrics_df$window == w, ]
    # best by each metric
    best_MAD <- m[order(m$MAD), ][1, ]
    best_IQR <- m[order(m$IQR), ][1, ]
    best_SD  <- m[order(m$SD),  ][1, ]
    best_CONS <- pick_consensus(m)

    best_rows <- rbind(best_rows,
      data.frame(fleet = fl, window = w, method = "MAD",  start_age = best_MAD$start_age,  end_age = best_MAD$end_age,  ages = best_MAD$ages,
                 MAD = best_MAD$MAD, IQR = best_MAD$IQR, SD = best_MAD$SD),
      data.frame(fleet = fl, window = w, method = "IQR",  start_age = best_IQR$start_age,  end_age = best_IQR$end_age,  ages = best_IQR$ages,
                 MAD = best_IQR$MAD, IQR = best_IQR$IQR, SD = best_IQR$SD),
      data.frame(fleet = fl, window = w, method = "SD",   start_age = best_SD$start_age,   end_age = best_SD$end_age,   ages = best_SD$ages,
                 MAD = best_SD$MAD, IQR = best_SD$IQR, SD = best_SD$SD),
      data.frame(fleet = fl, window = w, method = "CONSENSUS", start_age = best_CONS$start_age, end_age = best_CONS$end_age, ages = best_CONS$ages,
                 MAD = best_CONS$MAD, IQR = best_CONS$IQR, SD = best_CONS$SD)
    )
  }
}

write.csv(best_rows, file.path(outdir, "best_blocks_summary.csv"), row.names = FALSE)

# ---- SS-ready age bins (from consensus blocks) ----
# For SS, the age bins are simply the integer ages to include, e.g.:
#  13 14 15
# Here we report per fleet for each window.
ss_lines <- c("SS-ready age bins (3 consecutive ages) from CONSENSUS blocks",
              paste0("Input: ", infile),
              "")

for (fl in fleets) {
  ss_lines <- c(ss_lines, paste0("Fleet ", fl, ":"))
  for (w in windows) {
    row <- best_rows[best_rows$fleet == fl & best_rows$window == w & best_rows$method == "CONSENSUS", ][1, ]
    a <- row$start_age:row$end_age
    ss_lines <- c(ss_lines,
                  paste0("  window=", w, "  ages: ", paste(a, collapse = " ")))
  }
  ss_lines <- c(ss_lines, "")
}

writeLines(ss_lines, con = file.path(outdir, "SS_ready_age_bins.txt"))

# ---- plots: metric vs start_age ----
png_metric_plot <- function(m, fl, w, outdir) {
  fn <- file.path(outdir, paste0("fleet", fl, "_window", w, "_metrics.png"))
  png(fn, width = 1200, height = 800, res = 150)
  par(mfrow = c(2,2), mar = c(4,4,2,1))
  plot(m$start_age, m$MAD, type = "b", xlab = "Start age (block size = 3)", ylab = "MAD(residuals)",
       main = paste0("Fleet ", fl, " window=", w, " : MAD"))
  plot(m$start_age, m$IQR, type = "b", xlab = "Start age (block size = 3)", ylab = "IQR(residuals)",
       main = paste0("Fleet ", fl, " window=", w, " : IQR"))
  plot(m$start_age, m$SD, type = "b", xlab = "Start age (block size = 3)", ylab = "SD(residuals)",
       main = paste0("Fleet ", fl, " window=", w, " : SD"))
  plot(m$start_age, m$VAR, type = "b", xlab = "Start age (block size = 3)", ylab = "Var(residuals)",
       main = paste0("Fleet ", fl, " window=", w, " : VAR"))
  dev.off()
}

for (fl in fleets) {
  for (w in windows) {
    m <- metrics_df[metrics_df$fleet == fl & metrics_df$window == w, ]
    png_metric_plot(m, fl, w, outdir)
  }
}

cat("Done.\n")
cat("Wrote outputs to: ", normalizePath(outdir), "\n", sep = "")
