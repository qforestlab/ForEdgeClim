############################################################################################
# This script visualizes the parameter optimization process, using the CMA-ES algorithm (cmaesr).
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(ggh4x)
  library(grid)
})

start_analysis <- Sys.time()

# -----------------
# INPUT
# -----------------
season = "spring_50_" # winter spring summer autumn all_seasons

number_of_params = 3

# -----------------
# PLOTS — Convergence visuals
# -----------------

offspring_log <- read.csv(paste0("Output/calibration/data/logging_offspring_", season, "generations.csv"))
generation_log <- read.csv(paste0("Output/calibration/data/logging_generation_", season, "generations.csv"))
res <- readRDS(paste0("Output/calibration/data/calibration_results_", season, "generations.rds"))
output_path <- "Output/calibration/plots/"
lambda <- 4 + floor(3 * log(number_of_params))


# PLOT 1: Offspring RMSE between observed and modelled temperature
# ------

# Best offspring per generation + jitter-x
best_points <- offspring_log %>%
  group_by(iter) %>%
  slice_min(rmse, n = 1, with_ties = FALSE) %>%  # pick 1 value
  ungroup() %>%
  mutate(jitter_x = iter + (offspring_id - 1) / lambda)

p_RMSE <- ggplot() +
  # All offspring
  geom_point(data = offspring_log,
             aes(x = iter + (offspring_id - 1) / lambda, y = rmse, colour = "all offsprings")) +
  geom_smooth(data = offspring_log,
              aes(x = iter + (offspring_id - 1) / lambda, y = rmse, colour = "all offsprings"),
              se = FALSE, method = "loess", formula = y ~ x) +

  # Best offspring per generation is highlighted
  geom_point(data = best_points,
             aes(x = jitter_x, y = rmse, colour = "best offspring per generation"),
             size = 2) +
  geom_line(data = best_points,
            aes(x = jitter_x, y = rmse, colour = "best offspring per generation"),
            linewidth = 1, linetype = "dashed") +

  geom_vline(xintercept = unique(offspring_log$iter),
             linetype = "dashed", colour = "orange") +

  scale_colour_manual(values = c("all offsprings" = "blue",
                                 "best offspring per generation" = alpha("red", 0.5))) +
  guides(color = guide_legend(title = NULL, override.aes = list(size = 4, linewidth = 2))) +
  # Captions
  labs(x = "Generation",
       y = "RMSE (°C)",
       title = "(a) RMSE evolution over generations",
       colour = NULL) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 40),
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40)
  )

print(p_RMSE)
ggsave(paste0(output_path, "RMSE_", season, ".png"), plot = p_RMSE, width = 16, height = 10, dpi = 500)




# PLOT 2: 2D plot in PCA space
# ------

param_cols <- setdiff(names(offspring_log), c("eval_id","iter","offspring_id","rmse","time"))

# Data matrix (all generations together)
X <- offspring_log %>% select(all_of(param_cols)) %>% as.matrix()
keep <- complete.cases(X) # only keep complete parameter setes
Xc  <- scale(X[keep, ], center = TRUE, scale = TRUE)  # standardize/rescale parameters: mu = 0, variance = 1

iters_all <- offspring_log$iter[keep]

# PCA fit on all data (so axes stay consistent) (quite similar to PCs from CMA-ES algorithm: PCA on final population would be similar to final PCA in CMA-ES)
pca <- prcomp(Xc, center = FALSE, scale. = FALSE)
scores <- as.data.frame(pca$x[, 1:2]) # scores of the 2 main PCs
names(scores) <- c("PC1", "PC2")
scores$iter <- iters_all

# Facet plot per generation
p_PCA_space <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "blue") +
  stat_ellipse(type = "norm", level = 0.95, colour = "red", linewidth = 1) +
  coord_equal() +
  facet_wrap(~iter, ncol = 10) +
  labs(
    title = "(b) PCA on parameter space",
    #subtitle = "Points = offspring,\nEllipse = 95% spread",
    x = "PC1", y = "\nPC2"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 40),
    axis.title = element_text(size = 40),
    axis.text = element_text(size = 40),
    strip.text = element_text(size = 40) # facet labels
  ) +
  scale_x_continuous(labels = function(x) rep(" ", length(x))) +
  scale_y_continuous(labels = function(x) rep(" ", length(x)))


print(p_PCA_space)
ggsave(paste0(output_path, "PCA_space_", season, ".png"), plot = p_PCA_space, width = 16, height = 14, dpi = 500)



end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f s\n', as.numeric(end_analysis - start_analysis, units = 'secs')))


# Console summaries
cat("\nOptimised parameters (final):\n")
print(res$best.param)
cat("\nMinimal RMSE observed:\n")
print(res$best.fitness)
