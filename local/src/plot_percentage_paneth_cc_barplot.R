library(data.table)
library(mclust)
library(cluster)
library(reshape2)
library(ggplot2)
library(ramify)
library(psych)
library(pheatmap)
library(resample)
library(dplyr)
library(viridis)
library(igraph)
library(patchwork)
library(wesanderson)
library(reshape)
library(tidyr)

cellpress_theme <- function(base_size = 7) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.3),
      axis.ticks.length = unit(1.2, "mm"),
      axis.text = element_text(size = base_size),
      axis.title = element_text(size = base_size),
      plot.title = element_text(size = base_size, hjust = 0.5, face = "plain"),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size),
      plot.margin = margin(1.5, 1.5, 1.5, 1.5, unit = "mm")
    )
}

clust <- snakemake@input[["clust"]]
cc <- snakemake@input[["cc"]]
plot_out <- snakemake@output[["plot_out"]]
log_out <- snakemake@log[["log"]]

dato <- read.table(file = clust, row.names = 1, sep = ",", header = TRUE)

dato$isPaneth <- as.character(dato$isPaneth)

dato$isPaneth[dato$isPaneth == "filtered"] <- "Others"
dato$isPaneth[dato$isPaneth == "nPaneth"] <- "Others"

dato$isPaneth[dato$isPaneth == "Paneth"] <- "SPC"
dato$isPaneth[dato$isPaneth == "Others"] <- "NSPC"

dato$isPaneth <- factor(dato$isPaneth, levels = c("NSPC", "SPC"))

dato_cc <- read.table(file = cc, row.names = 1, sep = ",", header = FALSE)

dato <- merge(dato, dato_cc, by = "row.names")

dato$isPaneth <- factor(dato$isPaneth, levels = c("NSPC", "SPC"))
dato$V2 <- factor(dato$V2, levels = c("G1", "S", "G2M"))

dato.summary <- dato %>%
  group_by(isPaneth) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  mutate(percent = total_count / sum(total_count))

dato.summary_paneth <- dato %>%
  group_by(isPaneth, V2) %>%
  summarise(total_count = n(), .groups = "drop") %>%
  complete(
    isPaneth = factor(c("NSPC", "SPC"), levels = c("NSPC", "SPC")),
    V2 = factor(c("G1", "S", "G2M"), levels = c("G1", "S", "G2M")),
    fill = list(total_count = 0)
  ) %>%
  group_by(isPaneth) %>%
  mutate(
    percent = total_count / sum(total_count),
    mucca = 1 - (cumsum(percent) - 0.5 * percent)
  ) %>%
  ungroup()

colori_paneth <- c("NSPC" = "#76069A", "SPC" = "#099963")

df_paneth <- dato.summary %>%
  mutate(
    value = ifelse(isPaneth == "NSPC", -percent, percent),
    isPaneth = factor(isPaneth, levels = c("NSPC", "SPC"))
  )

paneth <- ggplot(df_paneth, aes(x = isPaneth, y = value, fill = isPaneth)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(
    aes(
      label = paste0(round(percent * 100, 1), "%"),
      hjust = ifelse(value > 0, -0.15, 1.15)
    ),
    size = 2.6
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5),
    labels = function(x) paste0(abs(x) * 100, "%"),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  scale_x_discrete(expand = expansion(add = 0)) +
  scale_fill_manual(values = colori_paneth) +
  labs(x = NULL, y = NULL) +
  cellpress_theme(base_size = 6) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    legend.spacing.y = unit(1, "mm"),
    legend.spacing.x = unit(1, "mm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank()
  )

# chi-square test: NSPC vs SPC across cell-cycle phases
cont <- dato.summary_paneth[, c("isPaneth", "V2", "total_count")] %>%
  mutate(
    isPaneth = factor(isPaneth, levels = c("NSPC", "SPC")),
    V2 = factor(V2, levels = c("G1", "S", "G2M"))
  )

tab <- xtabs(total_count ~ isPaneth + V2, data = cont)

chisq <- chisq.test(tab, correct = FALSE)
chisq_expected_low <- any(chisq$expected < 5)

df_plot <- cont %>%
  group_by(isPaneth) %>%
  mutate(percent = total_count / sum(total_count)) %>%
  ungroup() %>%
  mutate(
    value = ifelse(isPaneth == "NSPC", -percent, percent),
    V2 = factor(V2, levels = c("G1", "S", "G2M")),
    isPaneth = factor(isPaneth, levels = c("NSPC", "SPC"))
  )

colori_ciclo <- c("G1" = "#0173b4", "G2M" = "#93c2e9", "S" = "#b0d395")

ciclo <- ggplot(df_plot, aes(x = V2, y = value, fill = V2)) +
  geom_col(width = 0.9) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(
    aes(
      label = paste0(round(percent * 100, 1), "%"),
      hjust = ifelse(value > 0, -0.12, 1.12)
    ),
    size = 2.6
  ) +
  coord_flip(clip = "off") +
  scale_y_continuous(
    limits = c(-1, 1),
    breaks = seq(-1, 1, 0.5),
    labels = function(x) paste0(abs(x) * 100, "%"),
    expand = expansion(mult = c(0.02, 0.10))
  ) +
  scale_x_discrete(expand = expansion(add = 0)) +
  scale_fill_manual(values = colori_ciclo) +
  labs(x = NULL, y = "Percent of cells") +
  cellpress_theme(base_size = 6) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = c(0.98, 0.98),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 6),
    legend.key.size = unit(3, "mm"),
    legend.spacing.y = unit(1, "mm"),
    legend.spacing.x = unit(1, "mm"),
    legend.key = element_rect(fill = NA, colour = NA),
    legend.background = element_blank(),
    plot.margin = margin(1.5, 6, 1.5, 1.5, unit = "mm")
  )

l <- paneth / ciclo + plot_layout(heights = c(2, 3))

pdf(plot_out, width = 55 / 25.4, height = 55 / 25.4, useDingbats = FALSE)
print(l)
invisible(dev.off())

main <- function(log_f) {
  log_txt <- capture.output({

    cat("cellcycle percentage\n")
    print(dato.summary_paneth)

    cat("\n\nSPC / NSPC percentage\n")
    print(dato.summary)

    cat("\n\n==============================\n")
    cat("CHI-SQUARE TEST: CELL CYCLE DISTRIBUTION\n")
    cat("==============================\n")

    cat("\nTested hypothesis:\n")
    cat("H0: cell-cycle phase distribution is independent of SPC/NSPC identity.\n")
    cat("Groups: NSPC vs SPC\n")
    cat("Phases: G1, S, G2M\n")

    cat("\nObserved contingency table:\n")
    print(tab)

    cat("\nChi-square result:\n")
    cat(sprintf(
      "X-squared = %.4f, df = %d, p-value = %.3e\n",
      as.numeric(chisq$statistic),
      as.numeric(chisq$parameter),
      chisq$p.value
    ))

    if (chisq_expected_low) {
      cat("\nWARNING:\n")
      cat("At least one expected count is < 5. Chi-square approximation may be unreliable.\n")
    }

    cat("\nInterpretation:\n")
    cat("A significant p-value indicates that cell-cycle phase distribution differs between NSPC and SPC cells.\n")
  })

  writeLines(log_txt, con = log_f)
}

main(log_out)




