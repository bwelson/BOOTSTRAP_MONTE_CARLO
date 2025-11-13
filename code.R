# ============================================================================
# BOOTSTRAP CONFIDENCE INTERVAL METHODS FOR SMALL SAMPLES
# A Comprehensive Monte Carlo Study for African Research Contexts
# Author: Bentum Welson & Asante Akosua Agnes (A research project)
# Institution: KNUST, Ghana
# Date: March 2025 - August 2025
# Version: 3.0 (Enhanced)
# ============================================================================

# SETUP ======================================================================
rm(list = ls())
set.seed(23)

# Required packages
required_packages <- c("boot", "parallel", "ggplot2", "dplyr", "tidyr", 
                       "knitr", "kableExtra", "gridExtra", "patchwork",
                       "pheatmap", "viridis", "MASS", "microbenchmark")

new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

lapply(required_packages, library, character.only = TRUE)

theme_set(theme_minimal(base_size = 12))

# SIMULATION PARAMETERS ======================================================
cat("\n", rep("=", 80), "\n", sep="")
cat("ENHANCED BOOTSTRAP SIMULATION STUDY\n")
cat(rep("=", 80), "\n\n", sep="")

# Sample sizes
sample_sizes <- c(20, 30, 40, 50, 100, 200)

# Distributions (expanded)
distributions <- c("normal", "chisq", "exponential", "t", 
                   "lognormal", "contaminated_normal")

# Estimands (expanded)
estimands <- c("mean", "median", "sd", "cv")

# Number of simulations
n_simulations <- 10000  # Use 1000 for testing, 10000 for final
n_bootstrap <- 1000

# Confidence level
conf_level <- 0.95
alpha <- 1 - conf_level

# True population parameters
true_mean <- 5
true_median <- 5
true_sd <- 2
true_cv <- true_sd / true_mean

# Distribution parameters
dist_params <- list(
  normal = list(mean = true_mean, sd = true_sd),
  chisq = list(df = 10, shift = true_mean - 10),
  exponential = list(rate = 1/true_mean),
  t = list(df = 5, shift = true_mean),
  lognormal = list(meanlog = log(true_mean) - 0.5*log(1 + (true_sd/true_mean)^2),
                   sdlog = sqrt(log(1 + (true_sd/true_mean)^2))),
  contaminated_normal = list(mean = true_mean, sd = true_sd, 
                             contamination = 0.1, contam_scale = 5)
)

# True values for each estimand by distribution
true_values <- list(
  normal = list(mean = true_mean, median = true_mean, sd = true_sd, cv = true_cv),
  chisq = list(mean = true_mean, median = NA, sd = sqrt(10), cv = NA),
  exponential = list(mean = true_mean, median = true_mean * log(2), 
                     sd = true_mean, cv = 1),
  t = list(mean = true_mean, median = true_mean, sd = sqrt(5/3), cv = NA),
  lognormal = list(mean = true_mean, median = exp(log(true_mean) - 0.5*log(1.16)), 
                   sd = NA, cv = NA),
  contaminated_normal = list(mean = true_mean, median = true_mean, 
                             sd = true_sd, cv = true_cv)
)

# BOOTSTRAP METHODS ==========================================================

percentile_ci <- function(boot_obj, conf = 0.95) {
  alpha <- 1 - conf
  t_star <- boot_obj$t
  ci <- quantile(t_star, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)
  return(as.numeric(ci))
}

basic_ci <- function(boot_obj, conf = 0.95) {
  alpha <- 1 - conf
  t0 <- boot_obj$t0
  t_star <- boot_obj$t
  ci <- 2*t0 - quantile(t_star, probs = c(1 - alpha/2, alpha/2), na.rm = TRUE)
  return(as.numeric(ci))
}

normal_ci <- function(boot_obj, conf = 0.95) {
  t0 <- boot_obj$t0
  se <- sd(boot_obj$t, na.rm = TRUE)
  z <- qnorm(1 - (1-conf)/2)
  ci <- c(t0 - z*se, t0 + z*se)
  return(ci)
}

bca_ci <- function(boot_obj, conf = 0.95) {
  ci_obj <- tryCatch({
    boot.ci(boot_obj, conf = conf, type = "bca")
  }, error = function(e) NULL, warning = function(w) NULL)
  
  if(is.null(ci_obj) || is.null(ci_obj$bca)) {
    return(c(NA, NA))
  }
  return(ci_obj$bca[4:5])
}

# ASYMPTOTIC METHODS (for comparison) ========================================

t_interval <- function(data, conf = 0.95) {
  n <- length(data)
  m <- mean(data)
  se <- sd(data) / sqrt(n)
  t_crit <- qt(1 - (1-conf)/2, df = n-1)
  ci <- c(m - t_crit*se, m + t_crit*se)
  return(ci)
}

# STATISTIC FUNCTIONS ========================================================

mean_statistic <- function(data, indices) {
  return(mean(data[indices], na.rm = TRUE))
}

median_statistic <- function(data, indices) {
  return(median(data[indices], na.rm = TRUE))
}

sd_statistic <- function(data, indices) {
  return(sd(data[indices], na.rm = TRUE))
}

cv_statistic <- function(data, indices) {
  d <- data[indices]
  m <- mean(d, na.rm = TRUE)
  s <- sd(d, na.rm = TRUE)
  return(s / m)
}

get_statistic_function <- function(estimand) {
  switch(estimand,
         mean = mean_statistic,
         median = median_statistic,
         sd = sd_statistic,
         cv = cv_statistic,
         stop("Unknown estimand"))
}

get_true_value <- function(distribution, estimand) {
  true_values[[distribution]][[estimand]]
}

# DATA GENERATION ============================================================

generate_data <- function(n, distribution, params) {
  data <- switch(distribution,
    normal = rnorm(n, mean = params$mean, sd = params$sd),
    chisq = rchisq(n, df = params$df) + params$shift,
    exponential = rexp(n, rate = params$rate),
    t = rt(n, df = params$df) + params$shift,
    lognormal = rlnorm(n, meanlog = params$meanlog, sdlog = params$sdlog),
    contaminated_normal = {
      base <- rnorm(n, mean = params$mean, sd = params$sd)
      contam_idx <- rbinom(n, 1, params$contamination) == 1
      base[contam_idx] <- rnorm(sum(contam_idx), mean = params$mean, 
                                sd = params$sd * params$contam_scale)
      base
    },
    stop("Unknown distribution")
  )
  return(data)
}

# SINGLE REPLICATION =========================================================

run_single_replication <- function(n, distribution, params, estimand, 
                                   true_value, B = 1000) {
  # Generate data
  data <- generate_data(n, distribution, params)
  
  # Get statistic function
  stat_func <- get_statistic_function(estimand)
  
  # Run bootstrap
  boot_obj <- tryCatch({
    boot(data, statistic = stat_func, R = B)
  }, error = function(e) NULL)
  
  if(is.null(boot_obj)) return(NULL)
  
  # Calculate bootstrap CIs
  results <- tryCatch({
    list(
      estimand_value = boot_obj$t0,
      percentile = percentile_ci(boot_obj, conf = conf_level),
      basic = basic_ci(boot_obj, conf = conf_level),
      normal = normal_ci(boot_obj, conf = conf_level),
      bca = bca_ci(boot_obj, conf = conf_level),
      t_interval = if(estimand == "mean") t_interval(data, conf = conf_level) else c(NA, NA)
    )
  }, error = function(e) NULL)
  
  return(results)
}

# COVERAGE AND WIDTH =========================================================

check_coverage <- function(ci, true_value) {
  if(any(is.na(ci)) || is.na(true_value)) return(NA)
  return(ci[1] <= true_value & true_value <= ci[2])
}

calculate_width <- function(ci) {
  if(any(is.na(ci))) return(NA)
  return(ci[2] - ci[1])
}

check_tail_coverage <- function(ci, true_value) {
  if(any(is.na(ci)) || is.na(true_value)) return(list(left = NA, right = NA))
  list(
    left = ci[1] <= true_value,
    right = true_value <= ci[2]
  )
}

# MAIN SIMULATION ============================================================

run_simulation <- function(n, distribution, params, estimand, true_value,
                          n_sim, B, parallel = TRUE) {
  
  cat(sprintf("\nRunning: n=%d, dist=%s, estimand=%s, reps=%d\n", 
              n, distribution, estimand, n_sim))
  
  start_time <- Sys.time()
  
  if (parallel) {
    n_cores <- max(1, detectCores() - 1)
    cl <- makeCluster(n_cores)
    
    clusterExport(cl, varlist = c(
      "run_single_replication", "generate_data", 
      "mean_statistic", "median_statistic", "sd_statistic", "cv_statistic",
      "get_statistic_function", "percentile_ci", "basic_ci",
      "normal_ci", "bca_ci", "t_interval", "conf_level", "boot", 
      "n", "distribution", "params", "estimand", "true_value", "B"
    ), envir = environment())
    
    clusterEvalQ(cl, library(boot))
    
    results_list <- parLapply(cl, 1:n_sim, function(i) {
      run_single_replication(n, distribution, params, estimand, true_value, B)
    })
    
    stopCluster(cl)
    
  } else {
    results_list <- lapply(1:n_sim, function(i) {
      if (i %% 100 == 0) cat(sprintf("  Replication %d/%d\n", i, n_sim))
      run_single_replication(n, distribution, params, estimand, true_value, B)
    })
  }
  
  # Remove NULL results
  results_list <- results_list[!sapply(results_list, is.null)]
  n_successful <- length(results_list)
  
  cat(sprintf("  Successful: %d/%d (%.1f%%)\n", 
              n_successful, n_sim, 100 * n_successful / n_sim))
  
  # Process results
  methods <- c("percentile", "basic", "normal", "bca", "t_interval")
  summary_results <- data.frame(
    n = n,
    distribution = distribution,
    estimand = estimand,
    method = methods,
    coverage_rate = NA,
    mean_width = NA,
    sd_width = NA,
    left_tail_coverage = NA,
    right_tail_coverage = NA,
    n_successful = n_successful,
    stringsAsFactors = FALSE
  )
  
  for (method in methods) {
    cis <- t(sapply(results_list, function(x) x[[method]]))
    
    # Overall coverage
    coverage <- sapply(1:nrow(cis), function(i) 
      check_coverage(cis[i,], true_value))
    coverage_rate <- mean(coverage, na.rm = TRUE)
    
    # Tail coverage
    tail_cov <- lapply(1:nrow(cis), function(i) 
      check_tail_coverage(cis[i,], true_value))
    left_tail <- mean(sapply(tail_cov, function(x) x$left), na.rm = TRUE)
    right_tail <- mean(sapply(tail_cov, function(x) x$right), na.rm = TRUE)
    
    # Width statistics
    widths <- apply(cis, 1, calculate_width)
    mean_width <- mean(widths, na.rm = TRUE)
    sd_width <- sd(widths, na.rm = TRUE)
    
    # Store results
    idx <- summary_results$method == method
    summary_results[idx, "coverage_rate"] <- coverage_rate
    summary_results[idx, "mean_width"] <- mean_width
    summary_results[idx, "sd_width"] <- sd_width
    summary_results[idx, "left_tail_coverage"] <- left_tail
    summary_results[idx, "right_tail_coverage"] <- right_tail
  }
  
  # Add Monte Carlo standard errors
  summary_results$se_coverage <- sqrt(summary_results$coverage_rate * 
    (1 - summary_results$coverage_rate) / n_successful)
  
  elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("  Completed in %.1f minutes\n", elapsed))
  
  return(summary_results)
}

# RUN FULL SIMULATION STUDY ==================================================
cat("\n", rep("=", 80), "\n", sep="")
cat("STARTING FULL SIMULATION STUDY\n")
cat(rep("=", 80), "\n\n", sep="")

all_results <- data.frame()

for(dist in distributions) {
  for(estim in estimands) {
    # Skip if no true value
    true_val <- get_true_value(dist, estim)
    if(is.na(true_val)) next
    
    for(n in sample_sizes) {
      params <- dist_params[[dist]]
      
      result <- run_simulation(
        n = n,
        distribution = dist,
        params = params,
        estimand = estim,
        true_value = true_val,
        n_sim = n_simulations,
        B = n_bootstrap,
        parallel = TRUE
      )
      
      all_results <- rbind(all_results, result)
    }
  }
}

# Save results
save(all_results, file = "bootstrap_enhanced_results.RData")
write.csv(all_results, file = "bootstrap_enhanced_results.csv", row.names = FALSE)

cat("\n", rep("=", 80), "\n", sep="")
cat("SIMULATION COMPLETE!\n")
cat(rep("=", 80), "\n\n", sep="")

# COMPUTATIONAL TIMING ANALYSIS ==============================================
cat("\n=== COMPUTATIONAL TIMING ANALYSIS ===\n\n")

# Generate sample data for timing
timing_data <- rnorm(50, mean = 5, sd = 2)
timing_boot <- boot(timing_data, mean_statistic, R = 1000)

timing_results <- microbenchmark(
  percentile = percentile_ci(timing_boot),
  basic = basic_ci(timing_boot),
  normal = normal_ci(timing_boot),
  bca = bca_ci(timing_boot),
  times = 100
)

print(timing_results)
write.csv(summary(timing_results), "timing_analysis.csv", row.names = TRUE)

# REAL DATA ANALYSIS =========================================================
cat("\n=== REAL DATA EXAMPLES ===\n\n")

# Example 1: Small health study (simulated based on typical African data)
cat("Example 1: Maternal Health Study (n=25)\n")
maternal_health <- c(22.3, 24.1, 19.8, 25.6, 21.2, 23.4, 20.9, 26.1, 22.8, 24.5,
                     21.7, 23.9, 25.2, 20.4, 22.6, 24.8, 21.3, 23.1, 25.9, 22.4,
                     24.3, 21.8, 23.5, 20.6, 25.4)

real_boot1 <- boot(maternal_health, mean_statistic, R = 2000)
cat("Sample mean:", mean(maternal_health), "\n")
cat("Percentile CI:", percentile_ci(real_boot1), "\n")
cat("BCa CI:", bca_ci(real_boot1), "\n")
cat("t-interval:", t_interval(maternal_health), "\n\n")


# STATISTICAL TESTING ========================================================
cat("\n=== STATISTICAL SIGNIFICANCE TESTING ===\n\n")

# Test if methods differ significantly from 0.95
coverage_tests <- all_results %>%
  filter(estimand == "mean", n <= 50) %>%
  group_by(method) %>%
  summarise(
    mean_coverage = mean(coverage_rate, na.rm = TRUE),
    n_scenarios = n(),
    .groups = "drop"
  ) %>%
  mutate(
    # Binomial test for each method
    p_value = sapply(mean_coverage, function(cov) {
      binom.test(round(cov * 10000), 10000, p = 0.95)$p.value
    }),
    significant = p_value < 0.05
  )

print(coverage_tests)
write.csv(coverage_tests, "coverage_significance_tests.csv", row.names = FALSE)

# ENHANCED VISUALIZATIONS ====================================================
cat("\n=== CREATING ENHANCED VISUALIZATIONS ===\n\n")

# Focus on mean estimand for main plots
mean_results <- all_results %>% filter(estimand == "mean")

# 1. Coverage with Confidence Intervals
p1 <- ggplot(mean_results %>% filter(n <= 100), 
             aes(x = factor(n), y = coverage_rate, color = method)) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red", size = 1) +
  geom_point(size = 2.5, position = position_dodge(0.5)) +
  geom_errorbar(aes(ymin = coverage_rate - 1.96*se_coverage,
                    ymax = coverage_rate + 1.96*se_coverage),
                width = 0.3, position = position_dodge(0.5)) +
  facet_wrap(~distribution, nrow = 2) +
  labs(title = "Coverage Rates with 95% Monte Carlo Confidence Intervals",
       subtitle = "For mean estimator across distributions and sample sizes",
       x = "Sample Size", y = "Coverage Rate",
       color = "Method") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5))

library(tibble)
# 2. Heatmap of coverage rates
coverage_matrix <- mean_results %>%
  filter(method %in% c("percentile", "basic", "normal", "bca")) %>%
  mutate(scenario = paste0(distribution, "_n", n)) %>%
  dplyr::select(scenario, method, coverage_rate) %>%
  pivot_wider(names_from = method, values_from = coverage_rate) %>%
  column_to_rownames("scenario")


p2 <- pheatmap(coverage_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colorRampPalette(c("red", "yellow", "darkgreen"))(100),
         breaks = seq(0.88, 0.95, length.out = 101),
         main = "Coverage Rate Heatmap (Mean Estimator)",
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize_number = 8,
         silent = TRUE)

# 3. Tail Coverage Comparison
tail_data <- mean_results %>%
  filter(n <= 50, method %in% c("percentile", "basic", "normal", "bca")) %>%
  pivot_longer(cols = c(left_tail_coverage, right_tail_coverage),
               names_to = "tail", values_to = "coverage")

p3 <- ggplot(tail_data, aes(x = method, y = coverage, fill = tail)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.975, linetype = "dashed", color = "red") +
  facet_wrap(~paste0(distribution, " (n=", n, ")"), nrow = 3) +
  labs(title = "Left vs Right Tail Coverage (Small Samples)",
       subtitle = "Target: 0.975 for each tail",
       x = "Method", y = "Tail Coverage",
       fill = "Tail") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 4. Multi-estimand comparison
multi_est <- all_results %>%
  filter(n %in% c(30, 50), distribution == "normal",
         method %in% c("normal", "bca", "t_interval"))

p4 <- ggplot(multi_est, aes(x = estimand, y = coverage_rate, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "red") +
  facet_wrap(~paste0("n = ", n)) +
  labs(title = "Coverage Across Different Estimands (Normal Distribution)",
       x = "Estimand", y = "Coverage Rate",
       fill = "Method") +
  theme_minimal(base_size = 11)

# Save plots
ggsave("Enhanced_Fig1_Coverage_CI.png", plot = p1, 
       width = 14, height = 10, dpi = 300)
png("Enhanced_Fig2_Heatmap.png", width = 10, height = 8, units = "in", res = 300)
print(p2)
dev.off()
ggsave("Enhanced_Fig3_Tail_Coverage.png", plot = p3, 
       width = 14, height = 10, dpi = 300)
ggsave("Enhanced_Fig4_Multi_Estimand.png", plot = p4, 
       width = 10, height = 6, dpi = 300)

cat("Enhanced figures saved!\n")

# PUBLICATION TABLES =========================================================
cat("\n=== CREATING PUBLICATION TABLES ===\n\n")

# Table 1: Small sample performance (n ≤ 50, mean only)
table1 <- mean_results %>%
  filter(n <= 50, method != "t_interval") %>%
  group_by(method, distribution) %>%
  summarise(
    mean_cov = mean(coverage_rate, na.rm = TRUE),
    sd_cov = sd(coverage_rate, na.rm = TRUE),
    mean_width = mean(mean_width, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(distribution, desc(mean_cov))

write.csv(table1, "Enhanced_Table1_Small_Sample_Performance.csv", 
          row.names = FALSE)

# Table 2: Overall rankings
table2 <- mean_results %>%
  filter(method != "t_interval") %>%
  group_by(method) %>%
  summarise(
    scenarios = n(),
    mean_coverage = mean(coverage_rate, na.rm = TRUE),
    sd_coverage = sd(coverage_rate, na.rm = TRUE),
    coverage_below_94 = sum(coverage_rate < 0.94) / n(),
    mean_width = mean(mean_width, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_coverage))

write.csv(table2, "Enhanced_Table2_Overall_Rankings.csv", row.names = FALSE)

# Table 3: Recommendations by scenario
table3 <- mean_results %>%
  filter(n <= 50, method %in% c("normal", "bca")) %>%
  group_by(n, distribution) %>%
  summarise(
    best_method = method[which.max(coverage_rate)],
    best_coverage = max(coverage_rate),
    .groups = "drop"
  )

write.csv(table3, "Enhanced_Table3_Recommendations.csv", row.names = FALSE)

cat("Tables saved!\n")

# FINAL RECOMMENDATIONS ======================================================
cat("\n", rep("=", 80), "\n", sep="")
cat("PRACTICAL RECOMMENDATIONS FOR AFRICAN RESEARCHERS\n")
cat(rep("=", 80), "\n\n", sep="")

small_sample_performance <- mean_results %>%
  filter(n <= 50, method %in% c("percentile", "basic", "normal", "bca")) %>%
  group_by(method) %>%
  summarise(
    mean_coverage = mean(coverage_rate, na.rm = TRUE),
    coverage_close_to_95 = abs(mean_coverage - 0.95),
    mean_width = mean(mean_width, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(coverage_close_to_95, mean_width)

cat("For small samples (n ≤ 50):\n\n")
cat(sprintf("RECOMMENDED METHOD: %s\n", toupper(small_sample_performance$method[1])))
cat(sprintf("  - Coverage rate: %.4f (target: 0.9500)\n", 
            small_sample_performance$mean_coverage[1]))
cat(sprintf("  - Deviation from nominal: %.4f\n", 
            small_sample_performance$coverage_close_to_95[1]))
cat(sprintf("  - Mean CI width: %.4f\n\n", 
            small_sample_performance$mean_width[1]))

cat("ALTERNATIVE METHOD:", toupper(small_sample_performance$method[2]), "\n")
cat(sprintf("  - Coverage rate: %.4f\n", 
            small_sample_performance$mean_coverage[2]))
cat(sprintf("  - Mean CI width: %.4f\n\n", 
            small_sample_performance$mean_width[2]))

cat("\nGENERAL GUIDELINES:\n")
cat("1. Use BCa method when n ≥ 30 and computational time permits\n")
cat("2. Use Normal method when n ≥ 40 or for faster computation\n")
cat("3. Avoid Basic and Percentile methods for n < 50\n")
cat("4. Always use B ≥ 1000 bootstrap replications\n")
cat("5. Report both coverage performance and interval width\n")
cat("6. Consider classical t-intervals for comparison when n ≥ 30\n\n")

cat(rep("=", 80), "\n", sep="")
cat("STUDY COMPLETE! All results saved.\n")
cat(rep("=", 80), "\n", sep="")

