rm(list = ls())

library(tidyverse)

source("generate_randomness_fixed_covariates.R")
source("beta_optimization.R")

# set.seed(42)

omega_grid <- seq(0, 2, 0.1)


C1 <- 0.5
C2 <- 1
data <- generate_randomness_fixed_covariates(C1 = C1, C2 = C2)
tb_s1 <- map_dfr(
  omega_grid,
  function(omega){
    beta_optimization(
    data = data,
    gamma_C1 = 0.3, gamma_C2 = 1,
    gamma_A = -4,  
    omega = omega
    )
  }
) |>   
  mutate(
    subject = "S1",
    C1 = C1,
    C2 = C2
  )
tb_s1



C1 <- -0.5
C2 <- 1
data <- generate_randomness_fixed_covariates(C1 = C1, C2 = C2)
tb_s2 <- map_dfr(
  omega_grid,
  function(omega){
    beta_optimization(
    data = data,
    gamma_C1 = 0.3, gamma_C2 = 1,
    gamma_A = -4,
    omega = omega
    )
  }
) |>   
  mutate(
    subject = "S2",
    C1 = C1,
    C2 = C2
  )
tb_s2



C1 <- 0.5
C2 <- 0
data <- generate_randomness_fixed_covariates(C1 = C1, C2 = C2)
tb_s3 <- map_dfr(
  omega_grid,
  function(omega){
    beta_optimization(
    data = data,
    gamma_C1 = 0.3, gamma_C2 = 1,
    gamma_A = -4,  
    omega = omega
    )
  }
) |>   
  mutate(
    subject = "S3",
    C1 = C1,
    C2 = C2
  )
tb_s3



C1 <- -0.5
C2 <- 0
data <- generate_randomness_fixed_covariates(C1 = C1, C2 = C2)
tb_s4 <- map_dfr(
  omega_grid,
  function(omega){
    beta_optimization(
    data = data,
    gamma_C1 = 0.3, gamma_C2 = 1,
    gamma_A = -4,  
    omega = omega
    )
  }
) |>   
  mutate(
    subject = "S4",
    C1 = C1,
    C2 = C2
  )
tb_s4

figure_data <- bind_rows(tb_s1, tb_s2, tb_s3, tb_s4)

figure_data

p_beta <-
  figure_data |>
  ggplot(
    aes(
      x = omega,
      y = beta,
      color = subject,
      linetype = subject
    )
  ) +
  geom_line(linewidth = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = expression(omega),
    y = expression(beta)
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -6)
  )
p_beta
ggsave(
  "figures/fig_2_beta.pdf",
  plot = p_beta,
  width = 6,
  height = 4
)


p_Y <-
  figure_data |>
  ggplot(
    aes(
      x = omega,
      y = cost_Y,
      color = subject,
      linetype = subject
    )
  ) +
  geom_line(linewidth = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = expression(omega),
    y = "Biomarker risk"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -6)
  )
p_Y
ggsave(
  "figures/fig_2_Y.pdf",
  plot = p_Y,
  width = 6,
  height = 4
)

p_trt <-
  figure_data |>
  ggplot(
    aes(
      x = omega,
      y = cost_trt,
      color = subject,
      linetype = subject
    )
  ) +
  geom_line(linewidth = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = expression(omega),
    y = "Treatment risk"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -6)
  )
p_trt
ggsave(
  "figures/fig_2_trt.pdf",
  plot = p_trt,
  width = 6,
  height = 4
)


p_total <-
  figure_data |>
  ggplot(
    aes(
      x = omega,
      y = risk,
      color = subject,
      linetype = subject
    )
  ) +
  geom_line(linewidth = 2) +
  scale_color_brewer(palette = "Dark2") +
  labs(
    x = expression(omega),
    y = "Total risk"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -6)
  )
p_total
ggsave(
  "figures/fig_2_total.pdf",
  plot = p_total,
  width = 6,
  height = 4
)


p_pct <-
  figure_data |> 
  ggplot(
    aes(
      x = omega,
      y = pct_treated,
      color = subject,
      linetype = subject
    )
  ) +
  geom_line(linewidth = 2) +
  scale_color_brewer(palette = "Dark2") +
  # scale_y_continuous(labels = scales::percent_format()) +
  labs(
    x = expression(omega),
    y = "Proportion treated"
  ) +
  theme_minimal(base_size = 22) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -6)
  )
p_pct
ggsave(
  "figures/fig_2_pct.pdf",
  plot = p_pct,
  width = 6,
  height = 4
)
