rm(list = ls())

library(tidyverse)

source("generate_randomness_fixed_covariates.R")
source("risk_computation.R")

beta_grid <- seq(-40, 40, by = 0.25)
omegas <- c(0, 0.4, 0.6, 1, 3)


# Subject C1 = -0.5, C2 = 0
data <- generate_randomness_fixed_covariates(
  C1 = -0.5,
  C2 = 0
)

figure_data <- map_dfr(
  omegas,
  function(omega){
    map_dfr(
      beta_grid,
      function(beta){
        r <- risk_computation(          
          data = data,
          gamma_C1 = 0.3, gamma_C2 = 1,
          gamma_A = -4,  
          beta = beta,
          omega = omega
        )
        tibble(
          omega = omega,
          beta = beta,
          risk_Y = attr(r,"risk_Y"),
          risk_trt = attr(r,"risk_trt"),
          risk = r,
          pct_treated = attr(r,"pct_trt")
        )
      }
    )
  }
)

p_Y <- figure_data |>
  ggplot(aes(beta, risk_Y, group = factor(omega))) +
  geom_line(
    aes(
      linetype = factor(omega),
      color = factor(omega)
    ),
    linewidth = 2, 
    lineend = "round"
  ) + 
  scale_color_brewer(
    palette = "Dark2",
    name = expression(omega)
  ) +
  scale_linetype_discrete(
    name = expression(omega)
  ) + 
  labs(
    x = expression(beta),
    y = "Biomarker risk"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(),
    legend.margin = margin(t = -6)
  )
p_Y
ggsave(
  "figures/fig_2_sub_C1_m05_C2_0_Y.pdf",
  plot = p_Y,
  width = 6,
  height = 4
)

p_trt <- figure_data |>
  ggplot(aes(beta, risk_trt, group = factor(omega))) +
    geom_line(
      aes(
        linetype = factor(omega),
        color = factor(omega)
      ),
      linewidth = 2, 
      lineend = "round"
    ) +
    scale_color_brewer(
      palette = "Dark2",
      name = expression(omega)
    ) +
    scale_linetype_discrete(
      name = expression(omega)
    ) + 
    labs(
      x = expression(beta),
      y = "Treatment risk"
    ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(),
    legend.margin = margin(t = -6)
  )
p_trt
ggsave(
  "figures/fig_2_sub_C1_m05_C2_0_trt.pdf",
  plot = p_trt,
  width = 6,
  height = 4
)

p_risk <- figure_data |>
  ggplot(aes(beta, risk, group = factor(omega))) +
    geom_line(
      aes(
        linetype = factor(omega),
        color = factor(omega)
      ),
      linewidth = 2, 
      lineend = "round"
    ) +
    scale_color_brewer(
      palette = "Dark2",
      name = expression(omega)
    ) +
    scale_linetype_discrete(
      name = expression(omega)
    ) +
    labs(
      x = expression(beta),
      y = "Total risk",
      color = expression(omega),
      linetype = expression(omega)
    ) +
    theme_minimal(base_size = 22) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(),
      legend.margin = margin(t = -6)
    )
p_risk
ggsave(
  "figures/fig_2_sub_C1_m05_C2_0_risk.pdf",
  plot = p_risk,
  width = 6,
  height = 4
)

p_pct <- figure_data |> 
  ggplot(aes(beta, pct_treated, group = factor(omega)))  +
    geom_line(
      aes(
        linetype = factor(omega),
        color = factor(omega)
      ),
      linewidth = 2, 
      lineend = "round"
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_color_brewer(
      palette = "Dark2",
      name = expression(omega)
    ) +
    scale_linetype_discrete(
      name = expression(omega)
    ) + 
    labs(
      x = expression(beta),
      y = "Proportion treated",
      color = expression(omega),
      linetype = expression(omega)
    ) +
    theme_minimal(base_size = 22) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(),
      legend.margin = margin(t = -6)
    )
p_pct
ggsave(
  "figures/fig_2_sub_C1_m05_C2_0_pct.pdf",
  plot = p_pct,
  width = 6,
  height = 4
)

