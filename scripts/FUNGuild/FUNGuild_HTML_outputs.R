library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(plotly)
library(htmlwidgets)
library(limma)

project <- "/project/166726054/Marie-Eve_Analysis"

funout <- file.path(
  project,
  "08_functional_analysis/BITS/FUNGuild"
)

matrixdir <- file.path(funout, "matrices")
finaldir  <- file.path(funout, "final_tables")

outdir <- file.path(
  project,
  "10_figures/BITS_FUNGuild_HTML"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

metadata_file <- file.path(
  project,
  "06_exports/BITS_final/Metadata.tsv"
)

# ============================================================
# Helpers
# ============================================================

save_html <- function(widget, filename) {
  htmlwidgets::saveWidget(
    widget,
    file.path(outdir, filename),
    selfcontained = TRUE,
    libdir = NULL
  )
}

treat_order <- c(
  "CTRLMin",
  "CTROrg",
  "PGPRMin",
  "PGPROrg"
)

cult_order <- c(
  "Blueray",
  "Duke",
  "Patriot"
)

group_order <- as.vector(
  outer(
    treat_order,
    cult_order,
    paste,
    sep="__"
  )
)

# ============================================================
# 1. TROPHIC MODE HEATMAP
# ============================================================

tm <- read_tsv(
  file.path(
    finaldir,
    "all_assigned_trophic_mode_Treatment_x_Cultivar_mean.tsv"
  ),
  show_col_types=FALSE
)

tm_mat <- as.matrix(
  tm[, group_order]
)

rownames(tm_mat) <- tm$Category

# Convert to %
tm_pct <- 100 * tm_mat

heat_raw <- plot_ly(
  x = colnames(tm_pct),
  y = rownames(tm_pct),
  z = tm_pct,
  type = "heatmap",
  colorscale = "Viridis",
  hovertemplate =
    paste(
      "<b>%{y}</b>",
      "<br>Group: %{x}",
      "<br>% total fungal reads: %{z:.2f}",
      "<extra></extra>"
    )
) %>%
  layout(
    title = list(
      text = "FUNGuild trophic modes - Treatment x Cultivar"
    ),
    xaxis = list(
      title = "",
      tickangle = -45
    ),
    yaxis = list(
      title = "",
      automargin = TRUE
    ),
    margin = list(
      l = 270,
      b = 180,
      t = 80,
      r = 80
    )
  )

save_html(
  heat_raw,
  "FUNGuild_TrophicMode_heatmap_Treat_x_Cult.html"
)

# ============================================================
# 2. TROPHIC MODE Z-SCORE HEATMAP
# ============================================================

z <- t(scale(t(tm_pct)))
z[!is.finite(z)] <- 0

heat_z <- plot_ly(
  x = colnames(z),
  y = rownames(z),
  z = z,
  type = "heatmap",
  colorscale = list(
    c(0,   "#2166AC"),
    c(0.5, "#F7F7F7"),
    c(1,   "#B2182B")
  ),
  zmid = 0,
  hovertemplate =
    paste(
      "<b>%{y}</b>",
      "<br>Group: %{x}",
      "<br>Z-score: %{z:.2f}",
      "<extra></extra>"
    )
) %>%
  layout(
    title = list(
      text = "FUNGuild trophic modes - row Z-score"
    ),
    xaxis = list(
      title = "",
      tickangle = -45
    ),
    yaxis = list(
      title = "",
      automargin = TRUE
    ),
    margin = list(
      l = 270,
      b = 180,
      t = 80,
      r = 80
    )
  )

save_html(
  heat_z,
  "FUNGuild_TrophicMode_heatmap_Treat_x_Cult_Zscore.html"
)

# ============================================================
# 3. TROPHIC MODE - FERTILISATION
# ============================================================

tm_long <- tm %>%
  pivot_longer(
    cols = -Category,
    names_to = "Group",
    values_to = "Relative"
  ) %>%
  separate(
    Group,
    into = c("Traitement", "Cultivar"),
    sep = "__"
  ) %>%
  mutate(
    Fertilisation = ifelse(
      Traitement %in% c("CTRLMin","PGPRMin"),
      "Mineral",
      "Organic"
    ),
    PGPR = ifelse(
      Traitement %in% c("CTRLMin","CTROrg"),
      "No PGPR",
      "PGPR"
    )
  )

fert <- tm_long %>%
  group_by(
    Fertilisation,
    Category
  ) %>%
  summarise(
    Mean = mean(Relative),
    .groups="drop"
  )

write_tsv(
  fert,
  file.path(
    outdir,
    "FUNGuild_TrophicMode_Fertilisation_mean.tsv"
  )
)

p_fert <- ggplot(
  fert,
  aes(
    x = Fertilisation,
    y = 100 * Mean,
    fill = Category,
    text = paste0(
      "Trophic mode: ", Category,
      "<br>Fertilisation: ", Fertilisation,
      "<br>% total fungal reads: ",
      sprintf("%.2f", 100*Mean)
    )
  )
) +
  geom_col(width=0.72) +
  labs(
    title = "Predicted fungal trophic modes by fertilisation",
    x = NULL,
    y = "Relative abundance (% of total fungal reads)",
    fill = "Trophic mode"
  ) +
  theme_bw(base_size=12)

fert_html <- ggplotly(
  p_fert,
  tooltip="text"
)

save_html(
  fert_html,
  "FUNGuild_TrophicMode_Fertilisation.html"
)

# ============================================================
# 4. DIVERGING FERTILISATION EFFECT
# Organic - Mineral
# ============================================================

fert_diff <- fert %>%
  pivot_wider(
    names_from = Fertilisation,
    values_from = Mean
  ) %>%
  mutate(
    Difference = 100 * (Organic - Mineral)
  ) %>%
  arrange(Difference)

p_div <- ggplot(
  fert_diff,
  aes(
    x = Difference,
    y = reorder(Category, Difference),
    text = paste0(
      "Trophic mode: ", Category,
      "<br>Organic - Mineral: ",
      sprintf("%.2f", Difference),
      " percentage points"
    )
  )
) +
  geom_col() +
  geom_vline(
    xintercept=0,
    linetype=2
  ) +
  labs(
    title = "Fungal trophic-mode shifts associated with fertilisation",
    x = "Organic - Mineral (percentage points)",
    y = NULL
  ) +
  theme_bw(base_size=12)

div_html <- ggplotly(
  p_div,
  tooltip="text"
)

save_html(
  div_html,
  "FUNGuild_TrophicMode_Fertilisation_diverging_barplot.html"
)

# ============================================================
# 5. GUILD >=90% HEATMAP
# ============================================================

guild <- read_tsv(
  file.path(
    finaldir,
    "all_assigned_guild90_Treatment_x_Cultivar_mean.tsv"
  ),
  show_col_types=FALSE
)

guild_mat <- as.matrix(
  guild[, group_order]
)

rownames(guild_mat) <- guild$Category
guild_pct <- 100 * guild_mat

guild_heat <- plot_ly(
  x = colnames(guild_pct),
  y = rownames(guild_pct),
  z = guild_pct,
  type = "heatmap",
  colorscale = "Viridis",
  hovertemplate =
    paste(
      "<b>%{y}</b>",
      "<br>Group: %{x}",
      "<br>% total fungal reads: %{z:.2f}",
      "<extra></extra>"
    )
) %>%
  layout(
    title = list(
      text = "FUNGuild dominant guilds (>=90% cumulative abundance)"
    ),
    xaxis = list(
      title = "",
      tickangle = -45
    ),
    yaxis = list(
      title = "",
      automargin = TRUE
    ),
    margin = list(
      l = 560,
      b = 180,
      t = 90,
      r = 80
    )
  )

save_html(
  guild_heat,
  "FUNGuild_Guild90_heatmap_Treat_x_Cult.html"
)

# ============================================================
# 6. SAMPLE-LEVEL GUILD DIFFERENTIAL ANALYSIS
#    Fertilisation effect
# ============================================================

guild_sample <- read_tsv(
  file.path(
    matrixdir,
    "all_assigned_guild_relative_total_fungal.tsv"
  ),
  show_col_types=FALSE
)

meta <- read_tsv(
  metadata_file,
  comment="#",
  show_col_types=FALSE
)

sample_cols <- intersect(
  names(guild_sample),
  as.character(meta$`sample-id`)
)

meta <- meta %>%
  filter(
    as.character(`sample-id`) %in% sample_cols
  ) %>%
  mutate(
    SampleID = as.character(`sample-id`),
    Fertilisation = ifelse(
      Traitement %in% c("CTRLMin","PGPRMin"),
      "Mineral",
      "Organic"
    ),
    PGPR = ifelse(
      Traitement %in% c("CTRLMin","CTROrg"),
      "No",
      "Yes"
    ),
    Fertilisation = factor(
      Fertilisation,
      levels=c("Mineral","Organic")
    ),
    PGPR = factor(
      PGPR,
      levels=c("No","Yes")
    ),
    Cultivar = factor(Cultivar),
    Bloc = factor(Bloc)
  )

# Reorder matrix to metadata
sample_cols <- meta$SampleID

X <- as.matrix(
  guild_sample[, sample_cols]
)

rownames(X) <- guild_sample[[1]]

# ------------------------------------------------------------
# Prevalence filter:
# retain guilds observed in >= 12 / 48 samples
# (= 25% prevalence, same principle used previously)
# ------------------------------------------------------------

prev <- rowSums(X > 0)

keep <- prev >= 12

Xf <- X[keep, , drop=FALSE]

# logit-like stabilized transformation of proportions
# with sample-specific small pseudocount derived from
# smallest positive value, not an arbitrary fixed constant

minpos <- min(
  Xf[Xf > 0],
  na.rm=TRUE
)

eps <- minpos / 2

Y <- log2(
  (Xf + eps) /
  (1 - Xf + eps)
)

design <- model.matrix(
  ~ Fertilisation + PGPR + Cultivar + Bloc,
  data = meta
)

fit <- lmFit(
  Y,
  design
)

fit <- eBayes(fit)

coef_name <- grep(
  "^Fertilisation",
  colnames(design),
  value=TRUE
)[1]

tt <- topTable(
  fit,
  coef = coef_name,
  number = Inf,
  sort.by = "none"
)

tt$Guild <- rownames(tt)
tt$Prevalence <- prev[tt$Guild]

vol <- tt %>%
  transmute(
    Guild,
    log2FC = logFC,
    PValue = P.Value,
    FDR = adj.P.Val,
    Prevalence,
    NegLog10FDR = -log10(pmax(FDR, 1e-300)),
    Significant = FDR < 0.05
  )

write_tsv(
  vol,
  file.path(
    outdir,
    "FUNGuild_Guild_Fertilisation_limma.tsv"
  )
)

p_vol <- ggplot(
  vol,
  aes(
    x = log2FC,
    y = NegLog10FDR,
    text = paste0(
      "Guild: ", Guild,
      "<br>Effect Organic vs Mineral: ",
      sprintf("%.3f", log2FC),
      "<br>P = ", signif(PValue, 3),
      "<br>FDR = ", signif(FDR, 3),
      "<br>Prevalence = ", Prevalence, "/48"
    )
  )
) +
  geom_point(
    aes(shape=Significant),
    alpha=0.80,
    size=2.8
  ) +
  geom_hline(
    yintercept=-log10(0.05),
    linetype=2
  ) +
  geom_vline(
    xintercept=0,
    linetype=2
  ) +
  labs(
    title = "FUNGuild guild-level response to fertilisation",
    subtitle = "Positive effect = higher under organic fertilisation",
    x = "Estimated effect (Organic vs Mineral)",
    y = "-log10(FDR)",
    shape = "FDR < 0.05"
  ) +
  theme_bw(base_size=12)

vol_html <- ggplotly(
  p_vol,
  tooltip="text"
)

save_html(
  vol_html,
  "FUNGuild_Guild_Fertilisation_volcano.html"
)

# ============================================================
# Final QC
# ============================================================

cat("\n===== HTML OUTPUTS COMPLETE =====\n")

cat(
  "Guilds before prevalence filter:",
  nrow(X),
  "\n"
)

cat(
  "Guilds retained >=12/48 samples:",
  nrow(Xf),
  "\n"
)

cat(
  "Smallest positive relative abundance:",
  minpos,
  "\n"
)

cat(
  "Pseudocount used:",
  eps,
  "\n"
)

cat(
  "FDR < 0.05:",
  sum(vol$FDR < 0.05, na.rm=TRUE),
  "/",
  nrow(vol),
  "\n"
)

files <- list.files(
  outdir,
  pattern="\\.html$",
  full.names=FALSE
)

cat("\n")
cat(paste(sort(files), collapse="\n"))
cat("\n")
