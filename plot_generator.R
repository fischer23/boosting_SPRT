# This file creates all plots of the paper
# "Improving the (approximate) sequential probability ratio test by avoiding overshoot"

library(ggplot2)
library(patchwork)
library(dplyr)

### Figure 1 (simple null vs. simple alternative)

lab <- c("Boosted", "Conservative Power-one SPRT")
col <- c("cornflowerblue", "limegreen")


# load the data
load("results/simple.rda")

mus <- results_df$idx

# Calculate difference between stopping times in percent
max_val <- max(results_df$mean_stop_sprt)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt) / mean_stop_sprt)

# Plot for the mean sample size
p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_stop_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt, colour = "2", shape = "2")) +
  geom_text(aes(y = (mean_stop_boosted + mean_stop_sprt) / 2, label = sprintf("%.1f%%", percent_diff)),
    size = 3, nudge_x = 0.04, nudge_y = 4
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 16),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Strength of the signal") +
  ylab("Number of samples till rejection") +
  scale_x_continuous(breaks = mus, limits = c(min(mus) - 0.1, max(mus) + 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 300, 30), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )



# Plot for the mean type I error
p2 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_type_I_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_type_I_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_type_I_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_type_I_sprt, colour = "2", shape = "2")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 16),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Strength of the signal") +
  ylab("Type I error probability") +
  scale_x_continuous(breaks = mus, limits = c(min(mus) - 0.1, max(mus) + 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.06, 0.01), limits = c(0, 0.06), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

combined <- p1 + p2 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

# save the plot
ggsave("results/Plot_simple.pdf", plot = combined, width = 12, height = 4.5)


### Figure 2 (predictable plugin)

# load the data
load("results/comp_alt.rda")

mus <- results_df$idx

# Calculate difference between stopping times in percent
max_val <- max(results_df$mean_stop_sprt)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt) / mean_stop_sprt)

# Plot for the mean sample size
p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_stop_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt, colour = "2", shape = "2")) +
  geom_text(aes(y = (mean_stop_boosted + mean_stop_sprt) / 2, label = sprintf("%.1f%%", percent_diff)),
    size = 3, nudge_x = 0.04, nudge_y = 6
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 16),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Strength of the signal") +
  ylab("Number of samples till rejection") +
  scale_x_continuous(breaks = mus, limits = c(min(mus) - 0.1, max(mus) + 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 300, 30), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# save the plot
ggsave("results/Plot_comp_alt.pdf", plot = p1, width = 10, height = 4)


### Figure 3 (Sampling without replacement)

# load the data
load("results/WoR.rda")

lab <- c("Boosted", "RiLACS (Waudby-Smith et al. (2021))")
col <- c("cornflowerblue", "limegreen")

# Calculate difference between stopping times in percent
max_val <- max(results_df$mean_stop_sprt)
min_val <- min(results_df$mean_stop_boosted)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt) / mean_stop_sprt)

# Plot for the mean sample size
p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_stop_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt, colour = "2", shape = "2")) +
  geom_text(aes(y = (mean_stop_boosted + mean_stop_sprt) / 2, label = sprintf("%.1f%%", percent_diff)),
    size = 3, nudge_x = 0.1, nudge_y = -20
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 16),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Total number of samples") +
  ylab("Number of samples till rejection") +
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x)),
    limits = c(800, 150000)
  ) +
  scale_y_continuous(breaks = seq(0, max_val + 50, 100), limits = c(min_val - 50, max_val + 50), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# save the plot
ggsave("results/Plot_WoR.pdf", plot = p1, width = 10, height = 4)


### Figure 4 (stop for futility)

lab <- c("Boosted", "Wald's approximate SPRT", "Wald's conservative SPRT")
col <- c("cornflowerblue", "brown2", "limegreen")

# load the data
load("results/futility.rda")

betas <- results_df$idx

# Calculate difference between stopping times in percent with approximate SPRT
max_val <- max(results_df$mean_stop_sprt)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt) / mean_stop_sprt)

# Plot for the mean sample size with approximate SPRT
p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_stop_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt, colour = "2", shape = "2")) +
  geom_line(aes(y = -1, colour = "3", linetype = "3")) +
  geom_point(aes(y = -1, colour = "3", shape = "3")) +
  geom_text(aes(y = mean_stop_boosted, label = sprintf("%.1f%%", percent_diff)),
    size = 4, nudge_x = 0.005, nudge_y = -4
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2], "3" = col[3]),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17, "3" = 16),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  xlab(expression(beta)) +
  ylab("Numb. of samples till decision") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.03), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# Calculate difference between stopping times in percent with conservative SPRT
max_val <- max(results_df$mean_stop_sprt_cons)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt_cons) / mean_stop_sprt_cons)

# Plot for the mean sample size
p2 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = -1, colour = "2", linetype = "3")) +
  geom_point(aes(y = -1, colour = "2", shape = "2")) +
  geom_line(aes(y = mean_stop_sprt_cons, colour = "3", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt_cons, colour = "3", shape = "3")) +
  geom_text(aes(y = mean_stop_boosted, label = sprintf("%.1f%%", percent_diff)),
    size = 4, nudge_x = 0.005, nudge_y = -4
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2], "3" = col[3]),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17, "3" = 16),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  xlab(expression(beta)) +
  ylab("Numb. of samples till decision") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.03), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# Calculate difference between type II error probability in percent
p3 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = 1 - power_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = 1 - power_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = 1 - power_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = 1 - power_sprt, colour = "2", shape = "2")) +
  geom_line(aes(y = 1 - power_sprt_cons, colour = "3", linetype = "3")) +
  geom_point(aes(y = 1 - power_sprt_cons, colour = "3", shape = "3")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2], "3" = col[3]),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17, "3" = 16),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  xlab(expression(beta)) +
  ylab("Type II error probability") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.01), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0, 0.31), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# Plot for the mean type I error
p4 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_type_I_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_type_I_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_type_I_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_type_I_sprt, colour = "2", shape = "2")) +
  geom_line(aes(y = mean_type_I_sprt_cons, colour = "3", linetype = "3")) +
  geom_point(aes(y = mean_type_I_sprt_cons, colour = "3", shape = "3")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2], "3" = col[3]),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17, "3" = 16),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  xlab(expression(beta)) +
  ylab("Type I error probability") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.01), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.06, 0.01), limits = c(0, 0.06), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )


combined <- p1 + p2 + p3 + p4 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

# save the plot
ggsave("results/Plot_futility.pdf", plot = combined, width = 12, height = 8)


### Figure S.1 (Confidence sequence)

# load the data
load("results/CS.rda")

lab <- c("Boosted", "Robbins")
col <- c("cornflowerblue", "limegreen")

n <- length(results_df$idx)

p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = CS, colour = "2", linetype = "3")) +
  geom_line(aes(y = boosted_CS, colour = "1", linetype = "3")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Type", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Sample size") +
  ylab("Lower confidence bound") +
  scale_x_continuous(breaks = seq(0, n, 10), limits = c(0, n), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 2, 0.2), limits = c(0, 2), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# save the plot
ggsave("results/Plot_CS.pdf", plot = p1, width = 8, height = 4)


### Figure S.2 (simple null vs. simple alternative with alpha=0.01)

lab <- c("Boosted", "Conservative Power-one SPRT")
col <- c("cornflowerblue", "limegreen")


# load the data
load("results/simple_alpha001.rda")

mus <- results_df$idx

# Calculate difference between stopping times in percent
max_val <- max(results_df$mean_stop_sprt)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt) / mean_stop_sprt)

# Plot for the mean sample size
p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_stop_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt, colour = "2", shape = "2")) +
  geom_text(aes(y = (mean_stop_boosted + mean_stop_sprt) / 2, label = sprintf("%.1f%%", percent_diff)),
    size = 3, nudge_x = 0.04, nudge_y = 4
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 16),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Strength of the signal") +
  ylab("Number of samples till rejection") +
  scale_x_continuous(breaks = mus, limits = c(min(mus) - 0.1, max(mus) + 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 300, 30), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )



# Plot for the mean type I error
p2 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_type_I_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_type_I_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_type_I_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_type_I_sprt, colour = "2", shape = "2")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 16),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Strength of the signal") +
  ylab("Type I error probability") +
  scale_x_continuous(breaks = mus, limits = c(min(mus) - 0.1, max(mus) + 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.012, 0.002), limits = c(0, 0.012), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

combined <- p1 + p2 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

# save the plot
ggsave("results/Plot_simple_alpha001.pdf", plot = combined, width = 12, height = 4.5)


### Figure S.3 (stop for futility with alpha=0.01)

lab <- c("Boosted", "Wald's approximate SPRT", "Wald's conservative SPRT")
col <- c("cornflowerblue", "brown2", "limegreen")

# load the data
load("results/futility_alpha001.rda")

betas <- results_df$idx

# Calculate difference between stopping times in percent with approximate SPRT
max_val <- max(results_df$mean_stop_sprt)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt) / mean_stop_sprt)

# Plot for the mean sample size with approximate SPRT
p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_stop_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt, colour = "2", shape = "2")) +
  geom_line(aes(y = -1, colour = "3", linetype = "3")) +
  geom_point(aes(y = -1, colour = "3", shape = "3")) +
  geom_text(aes(y = mean_stop_boosted, label = sprintf("%.1f%%", percent_diff)),
    size = 4, nudge_x = 0.005, nudge_y = -5.5
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2], "3" = col[3]),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17, "3" = 16),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  xlab(expression(beta)) +
  ylab("Numb. of samples till decision") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.03), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# Calculate difference between stopping times in percent with conservative SPRT
max_val <- max(results_df$mean_stop_sprt_cons)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt_cons) / mean_stop_sprt_cons)

# Plot for the mean sample size
p2 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = -1, colour = "2", linetype = "3")) +
  geom_point(aes(y = -1, colour = "2", shape = "2")) +
  geom_line(aes(y = mean_stop_sprt_cons, colour = "3", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt_cons, colour = "3", shape = "3")) +
  geom_text(aes(y = mean_stop_boosted, label = sprintf("%.1f%%", percent_diff)),
    size = 4, nudge_x = 0.005, nudge_y = -5.5
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2], "3" = col[3]),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17, "3" = 16),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  xlab(expression(beta)) +
  ylab("Numb. of samples till decision") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.03), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# Calculate difference between type II error probability in percent
p3 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = 1 - power_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = 1 - power_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = 1 - power_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = 1 - power_sprt, colour = "2", shape = "2")) +
  geom_line(aes(y = 1 - power_sprt_cons, colour = "3", linetype = "3")) +
  geom_point(aes(y = 1 - power_sprt_cons, colour = "3", shape = "3")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2], "3" = col[3]),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17, "3" = 16),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  xlab(expression(beta)) +
  ylab("Type II error probability") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.01), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0, 0.31), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# Plot for the mean type I error
p4 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_type_I_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_type_I_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_type_I_sprt, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_type_I_sprt, colour = "2", shape = "2")) +
  geom_line(aes(y = mean_type_I_sprt_cons, colour = "3", linetype = "3")) +
  geom_point(aes(y = mean_type_I_sprt_cons, colour = "3", shape = "3")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2], "3" = col[3]),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17, "3" = 16),
    labels = c("1" = lab[1], "2" = lab[2], "3" = lab[3])
  ) +
  xlab(expression(beta)) +
  ylab("Type I error probability") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.01), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.012, 0.002), limits = c(0, 0.012), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )


combined <- p1 + p2 + p3 + p4 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

# save the plot
ggsave("results/Plot_futility_alpha001.pdf", plot = combined, width = 12, height = 8)


### Figure S.4 (simple null vs. simple alternative with Siegmund's SPRT)

lab <- c("Boosted", "Siegmund's power-one SPRT")
col <- c("cornflowerblue", "purple")


# load the data
load("results/simple.rda")

mus <- results_df$idx

# Calculate difference between stopping times in percent
max_val <- max(results_df$mean_stop_sprt_ds)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt_ds) / mean_stop_sprt_ds)

# Plot for the mean sample size
p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_stop_sprt_ds, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt_ds, colour = "2", shape = "2")) +
  geom_text(aes(y = (mean_stop_boosted + mean_stop_sprt_ds) / 2, label = sprintf("%.1f%%", percent_diff)),
    size = 3, nudge_x = 0.04, nudge_y = 4
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 16),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Strength of the signal") +
  ylab("Number of samples till rejection") +
  scale_x_continuous(breaks = mus, limits = c(min(mus) - 0.1, max(mus) + 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 300, 30), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )



# Plot for the mean type I error
p2 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_type_I_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_type_I_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_type_I_sprt_ds, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_type_I_sprt_ds, colour = "2", shape = "2")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 16),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab("Strength of the signal") +
  ylab("Type I error probability") +
  scale_x_continuous(breaks = mus, limits = c(min(mus) - 0.1, max(mus) + 0.1), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.06, 0.01), limits = c(0, 0.06), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

combined <- p1 + p2 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

# save the plot
ggsave("results/Plot_simple_Siegmund.pdf", plot = combined, width = 12, height = 4.5)


### Figure S.5 (stop for futility with Siegmund's SPRT)

lab <- c("Boosted", "Siegmund's SPRT")
col <- c("cornflowerblue", "purple")

# load the data
load("results/futility.rda")

betas <- results_df$idx

# Calculate difference between stopping times in percent with conservative SPRT
max_val <- max(results_df$mean_stop_sprt_ds)
results_df <- results_df %>%
  mutate(percent_diff = 100 * (mean_stop_boosted - mean_stop_sprt_ds) / mean_stop_sprt_ds)

# Plot for the mean sample size
p1 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_stop_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_stop_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_stop_sprt_ds, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_stop_sprt_ds, colour = "2", shape = "2")) +
  geom_text(aes(y = mean_stop_boosted, label = sprintf("%.1f%%", percent_diff)),
    size = 4, nudge_x = 0.01, nudge_y = 3.5
  ) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab(expression(beta)) +
  ylab("Numb. of samples till decision") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.03), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, max_val + 10), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# Calculate difference between type II error probability in percent
p2 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = 1 - power_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = 1 - power_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = 1 - power_sprt_ds, colour = "2", linetype = "3")) +
  geom_point(aes(y = 1 - power_sprt_ds, colour = "2", shape = "2")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  geom_abline(slope = 1, intercept = 0, color = "grey", linetype = "dashed") +
  xlab(expression(beta)) +
  ylab("Type II error probability") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.01), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.5, 0.05), limits = c(0, 0.31), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )

# Plot for the mean type I error
p3 <- ggplot(results_df, aes(idx)) +
  geom_line(aes(y = mean_type_I_boosted, colour = "1", linetype = "3")) +
  geom_point(aes(y = mean_type_I_boosted, colour = "1", shape = "1")) +
  geom_line(aes(y = mean_type_I_sprt_ds, colour = "2", linetype = "3")) +
  geom_point(aes(y = mean_type_I_sprt_ds, colour = "2", shape = "2")) +
  scale_linetype_manual(guide = "none", values = c("3" = "solid", "4" = "dashed", "5" = "dotted")) +
  scale_colour_manual(
    name = "Test", values = c("1" = col[1], "2" = col[2]),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  scale_shape_manual(
    name = "Test", values = c("1" = 15, "2" = 17),
    labels = c("1" = lab[1], "2" = lab[2])
  ) +
  xlab(expression(beta)) +
  ylab("Type I error probability") +
  scale_x_continuous(breaks = betas, limits = c(0, max(betas) + 0.01), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 0.06, 0.01), limits = c(0, 0.06), expand = c(0, 0)) +
  theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    axis.title.y = element_text(size = 15), axis.title.x = element_text(size = 15), axis.text = element_text(size = 10),
    legend.text = element_text(size = 15), legend.title = element_text(size = 15)
  )


combined <- p1 + p2 + p3 & theme(legend.position = "bottom")
combined <- combined + plot_layout(guides = "collect")

# save the plot
ggsave("results/Plot_futility_Siegmund.pdf", plot = combined, width = 14, height = 4.5)
