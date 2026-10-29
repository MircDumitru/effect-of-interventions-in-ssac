plot_trajectories <- function(tb, eta = 1.7, size = 12) {
  tb |> 
  ggplot(aes(time-1, Z, group = id)) +
    geom_line(
      color = "grey75", 
      alpha = 0.6,
      linewidth = 1.5, 
      lineend = "round"
    ) +
    geom_hline(
      yintercept = eta, 
      color = "#B2182B", 
      linetype = "22",
      linewidth = 2, 
      lineend = "round"
    ) +
    scale_x_continuous(
      breaks = unique(tb$time - 1)
    ) + 
    labs(
      x = "Time since baseline",
      y = "Z"
    ) +
    theme_minimal() +
    theme(
      axis.title = element_text(size = size),
      axis.text = element_text(size = size)
    ) 
}
