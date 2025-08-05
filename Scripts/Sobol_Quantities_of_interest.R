###############################################################################
# In this script all 'Sobol_all_at_once.R' rds outputfiles for the several
# seasons and time points are merged together to get an overview of the
# Quantities of Interest (4 seasons x 3 time points x 3 metrics).
###############################################################################


library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)
library(forcats)

#########
# INPUT #
#########

input_path <- "Output/sensitivity_analysis/Sobol_QoI/"
output_plot <- "Output/sensitivity_analysis/Sobol_QoI/Sobol_QoI.png"
output_plot_focused <- "Output/sensitivity_analysis/Sobol_QoI/Sobol_QoI_focused.png"

# File list with labels
files_info <- tribble(
  ~file, ~label,
  # # winter files
  #   # night
  # "400samples_25parameters_01h_27012025_avT.rds",     "Wi_Ni_≈T",
  # "400samples_25parameters_01h_27012025_SDT.rds",     "Wi_Ni_σT",
  # "400samples_25parameters_01h_27012025_gradT.rds",   "Wi_Ni_∇T",
  #   # morning
  # "400samples_25parameters_07h_27012025_avT.rds",     "Wi_Mo_≈T",
  # "400samples_25parameters_07h_27012025_SDT.rds",     "Wi_Mo_σT",
  # "400samples_25parameters_07h_27012025_gradT.rds",   "Wi_Mo_∇T",
  #   # noon
  # "400samples_25parameters_11h_27012025_avT.rds",     "Wi_No_≈T",
  # "400samples_25parameters_11h_27012025_SDT.rds",     "Wi_No_σT",
  # "400samples_25parameters_11h_27012025_gradT.rds",   "Wi_No_∇T",
  # # spring
  #   # night
  # "400samples_25parameters_01h_23042024_avT.rds",     "Sp_Ni_≈T",
  # "400samples_25parameters_01h_23042024_SDT.rds",     "Sp_Ni_σT",
  # "400samples_25parameters_01h_23042024_gradT.rds",   "Sp_Ni_∇T",
  #   # morning
  # "400samples_25parameters_06h_23042024_avT.rds",     "Sp_Mo_≈T",
  # "400samples_25parameters_06h_23042024_SDT.rds",     "Sp_Mo_σT",
  # "400samples_25parameters_06h_23042024_gradT.rds",   "Sp_Mo_∇T",
  #   # noon
  # "400samples_25parameters_12h_23042024_avT.rds",     "Sp_No_≈T",
  # "400samples_25parameters_12h_23042024_SDT.rds",     "Sp_No_σT",
  # "400samples_25parameters_12h_23042024_gradT.rds",   "Sp_No_∇T",
  # summer
    # night
  "400samples_25parameters_01h_07072023_avT.rds",     "Su_Ni_≈T",
  "400samples_25parameters_01h_07072023_SDT.rds",     "Su_Ni_σT",
  "400samples_25parameters_01h_07072023_gradT.rds",   "Su_Ni_∇T",
    # morning
  "400samples_25parameters_05h_07072023_avT.rds",     "Su_Mo_≈T",
  "400samples_25parameters_05h_07072023_SDT.rds",     "Su_Mo_σT",
  "400samples_25parameters_05h_07072023_gradT.rds",   "Su_Mo_∇T",
    # noon
  "400samples_25parameters_12h_07072023_avT.rds",     "Su_No_≈T",
  "400samples_25parameters_12h_07072023_SDT.rds",     "Su_No_σT",
  "400samples_25parameters_12h_07072023_gradT.rds",   "Su_No_∇T"#,
  # # autumn files
  #   # night
  # "400samples_25parameters_01h_05102024_avT.rds",     "Au_Ni_≈T",
  # "400samples_25parameters_01h_05102024_SDT.rds",     "Au_Ni_σT",
  # "400samples_25parameters_01h_05102024_gradT.rds",   "Au_Ni_∇T",
  #   # morning
  # "400samples_25parameters_07h_05102024_avT.rds",     "Au_Mo_≈T",
  # "400samples_25parameters_07h_05102024_SDT.rds",     "Au_Mo_σT",
  # "400samples_25parameters_07h_05102024_gradT.rds",   "Au_Mo_∇T",
  #   # noon
  # "400samples_25parameters_11h_05102024_avT.rds",     "Au_No_≈T",
  # "400samples_25parameters_11h_05102024_SDT.rds",     "Au_No_σT",
  # "400samples_25parameters_11h_05102024_gradT.rds",   "Au_No_∇T"

)


# Parameter sequence
param_order <- c(
  "betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v", "Kd_h", "Kb_h", "omega_g_h",
  "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v", "Kd_lw_h", "omega_g_lw_h",
  "h", "g_macro", "infl_macro", "infl_soil", "infl_forest", "g_forest", "p_ground", "g_soil", "k_soil"
)

# Color palette
param_colors <- c(
  colorRampPalette(c("gold", "darkorange"))(9),
  colorRampPalette(c("mediumpurple", "darkviolet"))(7),
  colorRampPalette(c("darkseagreen", "darkgreen"))(9)
)
names(param_colors) <- param_order

# Focused parameter sequence
focus_params <- c("h", "g_macro", "infl_macro", "infl_soil", "infl_forest",
                  "g_forest", "p_ground", "g_soil", "k_soil")

# Focused color palette: contrast rich
focus_colors <- setNames(RColorBrewer::brewer.pal(n = 9, name = "Paired"), focus_params)

# Load, combine and normalise
sobol_df <- files_info %>%
  mutate(data = map(file, ~ readRDS(file.path(input_path, .x)))) %>%
  mutate(S = map(data, ~ .x$S$original)) %>%
  select(label, S) %>%
  unnest(S) %>%
  mutate(parameter = rep(param_order, times = nrow(.) / length(param_order))) %>%
  group_by(label) %>%
  mutate(norm_value = pmax(0, S) / sum(pmax(0, S))) %>%
  ungroup() %>%
  mutate(parameter = factor(parameter, levels = param_order),
         label = fct_rev(fct_inorder(label)))

# Plot
p = ggplot(sobol_df, aes(x = norm_value, y = label, fill = parameter)) +
  geom_col(width = 0.5) +
  labs(
    title = "Parameter contribution to QoI variance (first-order Sobol indices)",
    x = "Normalised first-order Sobol-index",
    y = "QoI",
    fill = "Parameter",
    caption = "QoI = Quantity of Interest\nSp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  scale_fill_manual(values = param_colors) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        plot.caption = element_text(hjust = 0)
        )
ggsave(output_plot, plot = p, width = 10, height = 10, dpi = 500)

# Focused output plot

# Filter sobol_df on focused parameters
sobol_focus_df <- sobol_df %>%
  filter(parameter %in% focus_params)

# Focused plot
p_focus <- ggplot(sobol_focus_df, aes(x = norm_value, y = label, fill = parameter)) +
  geom_col(width = 0.5) +
  labs(
    title = "Key parameters contribution to QoI variance (first-order Sobol indices)",
    x = "Normalised first-order Sobol-index",
    y = "QoI",
    fill = "Parameter",
    caption = "QoI = Quantity of Interest\nSp = spring | Su = summer | Au = autumn | Wi = winter\nMo = morning | No = noon | Ni = night\n≈T = average temperature | σT = standard deviation on temperature | ∇T = temperature gradient"
  ) +
  scale_fill_manual(values = focus_colors) +
  theme_bw(base_size = 14) +
  theme(plot.title = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        plot.caption = element_text(hjust = 0)
  )

ggsave(output_plot_focused, plot = p_focus, width = 10, height = 10, dpi = 500)


