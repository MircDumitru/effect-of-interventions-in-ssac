rm(list = ls())

library(tidyverse)


# Subject and model parameters
# -------------------------------------------------------------
mu1_i <- 0.5
C1_i <- 0.5
C2_i <- 1

gamma_C1 <- 0.3
gamma_C2 <- 1
gamma_A  <- -4
# -------------------------------------------------------------


# Fixed reference proportion of subject-visits under treatment
# -------------------------------------------------------------
p_treated <- 0.53

d0 <- mu1_i +
  gamma_C1 * C1_i +
  gamma_C2 * C2_i

d1 <- d0 + gamma_A

signal_power <- (1 - p_treated) * d0^2 +
  p_treated * d1^2

cat("Untreated drift:", d0, "\n")
cat("Treated drift:", d1, "\n")
cat("Reference signal power:", signal_power, "\n")
# -------------------------------------------------------------


# SNR function
# -------------------------------------------------------------
snr_dynamic_db <- function(
    tau,
    sigma_eps,
    signal_power
) {
  10 * log10(
    signal_power /
      (tau^2 + 2 * sigma_eps^2)
  )
}
# -------------------------------------------------------------


# Heat-map grid
# -------------------------------------------------------------
tau_grid <- seq(
  from = 0.1,
  to = 2.3,
  length.out = 500
)

sigma_eps_grid <- seq(
  from = 0.02,
  to = 0.60,
  length.out = 500
)

snr_heatmap_data <- expand_grid(
  tau = tau_grid,
  sigma_eps = sigma_eps_grid
) |>
  mutate(
    snr_db = snr_dynamic_db(
      tau = tau,
      sigma_eps = sigma_eps,
      signal_power = signal_power
    )
  )



# Three configurations used in the article
# -------------------------------------------------------------
article_noise_settings <- tribble(
  ~setting,    ~tau, ~sigma_eps,
  "Low",        0.6,       0.15,
  "Moderate",   1.0,       0.25,
  "High",       2.0,       0.50
) |>
  mutate(
    snr_db = snr_dynamic_db(
      tau = tau,
      sigma_eps = sigma_eps,
      signal_power = signal_power
    ),
    snr_label = sprintf("%.1f dB", snr_db)
  )

print(article_noise_settings)



# Five target-SNR configurations used in the MC study
# -------------------------------------------------------------
target_snr <- c(0, 5, 10, 15, 20)
tau_eps_ratio <- 4

mc_noise_settings <- tibble(
  target_snr = target_snr
) |>
  mutate(
    sigma_eps = sqrt(
      signal_power /
        (
          (tau_eps_ratio^2 + 2) *
            10^(target_snr / 10)
        )
    ),
    tau = tau_eps_ratio * sigma_eps,
    snr_label = paste0(target_snr, " dB")
  )

print(mc_noise_settings)



# Heat map
# -------------------------------------------------------------
p_snr_heatmap <- ggplot() +
  geom_raster(
    data = snr_heatmap_data,
    aes(
      x = tau,
      y = sigma_eps,
      fill = snr_db
    ),
    interpolate = FALSE
  ) +
  geom_contour(
    data = snr_heatmap_data,
    aes(
      x = tau,
      y = sigma_eps,
      z = snr_db
    ),
    breaks = c(0, 5, 10, 15, 20),
    colour = "white",
    linewidth = 0.55,
    alpha = 0.85
  ) +
  geom_point(
    data = article_noise_settings,
    aes(
      x = tau,
      y = sigma_eps
    ),
    inherit.aes = FALSE,
    shape = 21,
    fill = "white",
    colour = "black",
    stroke = 1.1,
    size = 4
  ) +
  geom_label(
    data = article_noise_settings,
    aes(
      x = tau,
      y = sigma_eps,
      # label = paste0(setting, "\n", snr_label)
      label = snr_label
    ),
    nudge_y = 0.05,
    inherit.aes = FALSE,
    size = 3.8,
    fontface = "bold",
    colour = "black",
    fill = scales::alpha("white", 0.85),
    linewidth = 0.25,
    label.padding = unit(0.15, "lines")
  ) +
  scale_fill_viridis_c(
    option = "D",
    direction = 1,
    name = "SNR (dB)"
  ) +
  scale_x_continuous(
    breaks = seq(0, 2.0, by = 0.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = seq(0, 0.6, by = 0.1),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    xlim = c(0.1, 2.3),
    ylim = c(0.02, 0.60),
    clip = "off"
  ) +
  labs(
    x = expression(tau),
    y = expression(sigma[epsilon])
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = grid::unit(2, "cm"),
    panel.grid = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )
p_snr_heatmap

# ggsave(
#   "figures/web_fig_1_heatmap_dynamic_snr.png",
#   p_snr_heatmap,
#   width = 6,
#   height = 4,
#   units = "in",
#   dpi = 600,
#   bg = "white"
# )

ggsave(
    "figures/web_fig_1_heatmap_dynamic_snr.pdf",
    p_snr_heatmap,
    width = 6,
    height = 4
)


