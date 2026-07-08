suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(RColorBrewer)
})

metadata <- read.delim(
  "combined_metadata_three_controls.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)
feature_table <- read.delim(
  "combined_metagenome_three_controls.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)

names(metadata) <- trimws(names(metadata))
metadata[] <- lapply(metadata, function(x) if (is.character(x)) trimws(x) else x)
names(feature_table)[1] <- "Feature"

sample_ids <- intersect(names(feature_table)[-1], metadata[["#NAME"]])
metadata <- metadata[match(sample_ids, metadata[["#NAME"]]), , drop = FALSE]

abundance <- t(as.matrix(feature_table[, sample_ids, drop = FALSE]))
storage.mode(abundance) <- "numeric"

dir.create("S2a.analysis_outputs", showWarnings = FALSE)

make_plot <- function(cohort_name, healthy_n_expected) {
  cohort_meta <- metadata[metadata$Cohort == cohort_name, , drop = FALSE]

  cohort_meta$PCoA_group <- NA_character_
  cohort_meta$PCoA_group[cohort_meta$Group == "Healthy"] <- "Healthy"
  cohort_meta$PCoA_group[cohort_meta$new_grouping_info_cutoff10 == "Pa"] <- "PAD"
  cohort_meta$PCoA_group[cohort_meta$new_grouping_info_cutoff10 == "Hi"] <- "HID"
  cohort_meta$PCoA_group[cohort_meta$new_grouping_info_cutoff10 == "PPM"] <- "PPMD"
  cohort_meta$PCoA_group[cohort_meta$new_grouping_info_cutoff10 == "Commensal"] <- "NPD"

  group_levels <- c("PAD", "HID", "PPMD", "NPD", "Healthy")
  cohort_meta <- cohort_meta[
    !is.na(cohort_meta$PCoA_group) & cohort_meta$PCoA_group %in% group_levels,
    ,
    drop = FALSE
  ]
  cohort_meta$PCoA_group <- factor(cohort_meta$PCoA_group, levels = group_levels)

  cohort_abundance <- abundance[cohort_meta[["#NAME"]], , drop = FALSE]
  bray_dist <- vegdist(cohort_abundance, method = "bray")
  pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2)
  variance <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 2)

  plot_data <- data.frame(
    SampleID = cohort_meta[["#NAME"]],
    Group = cohort_meta$PCoA_group,
    PCoA1 = pcoa$points[, 1],
    PCoA2 = pcoa$points[, 2],
    stringsAsFactors = FALSE
  )

  adonis_fit <- adonis2(bray_dist ~ PCoA_group, data = cohort_meta, permutations = 999)
  adonis_r2 <- adonis_fit$R2[1]
  adonis_p <- adonis_fit$`Pr(>F)`[1]
  adonis_label <- sprintf("Adonis R2 = %.3f\nP = %.3g", adonis_r2, adonis_p)

  sample_counts <- as.data.frame(table(plot_data$Group), stringsAsFactors = FALSE)
  names(sample_counts) <- c("Group", "n")
  write.table(
    sample_counts,
    file = sprintf("S2a.analysis_outputs/%s_PCoA_sample_counts.tsv", cohort_name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  if (sample_counts$n[sample_counts$Group == "Healthy"] != healthy_n_expected) {
    warning(sprintf(
      "%s Healthy sample count is %s, expected %s.",
      cohort_name,
      sample_counts$n[sample_counts$Group == "Healthy"],
      healthy_n_expected
    ))
  }

  set1_colors <- brewer.pal(5, "Set1")
  names(set1_colors) <- group_levels

  p <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2, color = Group)) +
    geom_point(size = 2.8, alpha = 0.85) +
    stat_ellipse(aes(group = Group), linewidth = 0.5, linetype = 2, show.legend = FALSE) +
    scale_color_manual(values = set1_colors, drop = FALSE) +
    labs(
      title = sprintf("%s PCoA", cohort_name),
      x = sprintf("PCoA1 (%.2f%%)", variance[1]),
      y = sprintf("PCoA2 (%.2f%%)", variance[2]),
      color = NULL
    ) +
    annotate(
      "text",
      x = Inf,
      y = Inf,
      label = adonis_label,
      hjust = 1.05,
      vjust = 1.2,
      size = 4
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    )

  ggsave(
    filename = sprintf("S2a.analysis_outputs/%s_PCoA_Set1_Adonis.pdf", cohort_name),
    plot = p,
    width = 6,
    height = 5
  )
  ggsave(
    filename = sprintf("S2a.analysis_outputs/%s_PCoA_Set1_Adonis.png", cohort_name),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )

  adonis_out <- data.frame(
    Cohort = cohort_name,
    Term = rownames(adonis_fit),
    Df = adonis_fit$Df,
    SumOfSqs = adonis_fit$SumOfSqs,
    R2 = adonis_fit$R2,
    F = adonis_fit$F,
    P = adonis_fit$`Pr(>F)`,
    row.names = NULL,
    check.names = FALSE
  )
  write.table(
    adonis_out,
    file = sprintf("S2a.analysis_outputs/%s_PCoA_adonis.tsv", cohort_name),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

make_plot("CAMEB2", healthy_n_expected = 25)
make_plot("EMBARC", healthy_n_expected = 88)

message("Done. Wrote PCoA plots and Adonis results to S2a.analysis_outputs/")
