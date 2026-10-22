library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(pheatmap)

project <- "/project/166726054/Marie-Eve_Analysis"

infile <- file.path(
  project,
  "08_functional_analysis/BITS/FUNGuild/final_tables",
  "all_assigned_trophic_mode_Treatment_x_Cultivar_mean.tsv"
)

figdir <- file.path(
  project,
  "10_figures/BITS_FUNGuild"
)

dir.create(figdir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Read table
# ------------------------------------------------------------
dat <- read_tsv(
  infile,
  show_col_types = FALSE
)

mat <- as.data.frame(dat)

rownames(mat) <- mat$Category
mat$Category <- NULL

mat <- as.matrix(mat)

# Explicit experimental order
treatment_order <- c(
  "CTRLMin",
  "CTROrg",
  "PGPRMin",
  "PGPROrg"
)

cultivar_order <- c(
  "Blueray",
  "Duke",
  "Patriot"
)

desired_columns <- as.vector(
  outer(
    treatment_order,
    cultivar_order,
    paste,
    sep = "__"
  )
)

mat <- mat[, desired_columns, drop = FALSE]

cat("===== INPUT MATRIX =====\n")
cat("Rows       :", nrow(mat), "\n")
cat("Columns    :", ncol(mat), "\n")
cat("Min value  :", min(mat), "\n")
cat("Max value  :", max(mat), "\n")

# ------------------------------------------------------------
# 2. Column annotation
# ------------------------------------------------------------
annotation_col <- data.frame(
  Treatment = sub("__.*$", "", colnames(mat)),
  Cultivar  = sub("^.*__", "", colnames(mat))
)

rownames(annotation_col) <- colnames(mat)

# ------------------------------------------------------------
# 3. RAW HEATMAP
# Values = fraction of total fungal reads
# ------------------------------------------------------------
pdf(
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Treat_x_Cult_heatmap.pdf"
  ),
  width = 12,
  height = 7
)

pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 8,
  main = "FUNGuild trophic modes – relative abundance"
)

dev.off()

png(
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Treat_x_Cult_heatmap.png"
  ),
  width = 3600,
  height = 2100,
  res = 300
)

pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 8,
  main = "FUNGuild trophic modes – relative abundance"
)

dev.off()

# ------------------------------------------------------------
# 4. ROW Z-SCORE HEATMAP
# ------------------------------------------------------------
zmat <- t(scale(t(mat)))

# Replace undefined Z-scores if a row has zero variance
zmat[!is.finite(zmat)] <- 0

pdf(
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Treat_x_Cult_Zscore_heatmap.pdf"
  ),
  width = 12,
  height = 7
)

pheatmap(
  zmat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 8,
  main = "FUNGuild trophic modes – row Z-score"
)

dev.off()

png(
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Treat_x_Cult_Zscore_heatmap.png"
  ),
  width = 3600,
  height = 2100,
  res = 300
)

pheatmap(
  zmat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  border_color = NA,
  fontsize_row = 9,
  fontsize_col = 8,
  main = "FUNGuild trophic modes – row Z-score"
)

dev.off()

# ------------------------------------------------------------
# 5. Long-format table
# ------------------------------------------------------------
long <- dat %>%
  pivot_longer(
    cols = -Category,
    names_to = "Group",
    values_to = "RelativeAbundance"
  ) %>%
  separate(
    Group,
    into = c("Traitement", "Cultivar"),
    sep = "__",
    remove = FALSE
  )

long$Traitement <- factor(
  long$Traitement,
  levels = treatment_order
)

long$Cultivar <- factor(
  long$Cultivar,
  levels = cultivar_order
)

# ------------------------------------------------------------
# 6. Mean by treatment
# Equal design: 3 cultivars × 4 blocks
# ------------------------------------------------------------
treatment_summary <- long %>%
  group_by(Traitement, Category) %>%
  summarise(
    MeanRelativeAbundance = mean(RelativeAbundance),
    .groups = "drop"
  )

write_tsv(
  treatment_summary,
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Treatment_mean.tsv"
  )
)

p1 <- ggplot(
  treatment_summary,
  aes(
    x = Traitement,
    y = MeanRelativeAbundance,
    fill = Category
  )
) +
  geom_col(width = 0.75) +
  labs(
    x = "Treatment",
    y = "Relative abundance\n(fraction of total fungal reads)",
    fill = "Trophic mode",
    title = "Predicted fungal trophic modes by treatment"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    legend.position = "right"
  )

ggsave(
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Treatment_stacked_barplot.pdf"
  ),
  p1,
  width = 10,
  height = 7
)

ggsave(
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Treatment_stacked_barplot.png"
  ),
  p1,
  width = 10,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 7. Fertilisation × PGPR
# ------------------------------------------------------------
fp <- long %>%
  mutate(
    Fertilisation = case_when(
      Traitement %in% c("CTRLMin", "PGPRMin") ~ "Mineral",
      Traitement %in% c("CTROrg", "PGPROrg") ~ "Organic"
    ),
    PGPR = case_when(
      Traitement %in% c("CTRLMin", "CTROrg") ~ "No PGPR",
      Traitement %in% c("PGPRMin", "PGPROrg") ~ "PGPR"
    )
  ) %>%
  group_by(
    Fertilisation,
    PGPR,
    Category
  ) %>%
  summarise(
    MeanRelativeAbundance = mean(RelativeAbundance),
    .groups = "drop"
  )

write_tsv(
  fp,
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Fertilisation_PGPR_mean.tsv"
  )
)

fp$Fertilisation <- factor(
  fp$Fertilisation,
  levels = c("Mineral", "Organic")
)

fp$PGPR <- factor(
  fp$PGPR,
  levels = c("No PGPR", "PGPR")
)

p2 <- ggplot(
  fp,
  aes(
    x = PGPR,
    y = MeanRelativeAbundance,
    fill = Category
  )
) +
  geom_col(width = 0.75) +
  facet_wrap(~ Fertilisation) +
  labs(
    x = NULL,
    y = "Relative abundance\n(fraction of total fungal reads)",
    fill = "Trophic mode",
    title = "Predicted fungal trophic modes: fertilisation × PGPR"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right"
  )

ggsave(
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Fertilisation_PGPR_barplot.pdf"
  ),
  p2,
  width = 11,
  height = 7
)

ggsave(
  file.path(
    figdir,
    "FUNGuild_TrophicMode_Fertilisation_PGPR_barplot.png"
  ),
  p2,
  width = 11,
  height = 7,
  dpi = 300
)

# ------------------------------------------------------------
# 8. Final report
# ------------------------------------------------------------
cat("\n===== FIGURES COMPLETE =====\n")

files <- list.files(
  figdir,
  pattern = "FUNGuild_TrophicMode",
  full.names = FALSE
)

cat(paste(sort(files), collapse = "\n"))
cat("\n")
