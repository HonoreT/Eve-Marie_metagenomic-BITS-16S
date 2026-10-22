library(readr)
library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)

project <- "/project/166726054/Marie-Eve_Analysis"

matrix_file <- file.path(
  project,
  "08_functional_analysis/EMP/targeted_modules",
  "EMP_core_functions_relative.tsv"
)

stats_file <- file.path(
  project,
  "08_functional_analysis/EMP/targeted_modules",
  "EMP_core_functions_limma_ALL.tsv"
)

metadata_file <- file.path(
  project,
  "06_exports/EMP_final/Metadata.tsv"
)

outdir <- file.path(
  project,
  "10_figures/EMP_targeted_HTML"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. LOAD
# ============================================================

X <- read_tsv(
  matrix_file,
  show_col_types = FALSE
)

stats <- read_tsv(
  stats_file,
  show_col_types = FALSE
)

meta <- read_tsv(
  metadata_file,
  comment = "#",
  show_col_types = FALSE
)

meta <- meta %>%
  mutate(
    SampleID = as.character(`sample-id`),

    Fertilisation = case_when(
      Traitement %in% c("CTRLMin", "PGPRMin") ~ "Mineral",
      Traitement %in% c("CTROrg", "PGPROrg") ~ "Organic"
    ),

    PGPR = case_when(
      Traitement %in% c("CTRLMin", "CTROrg") ~ "No",
      Traitement %in% c("PGPRMin", "PGPROrg") ~ "Yes"
    ),

    Fertilisation = factor(
      Fertilisation,
      levels = c("Mineral", "Organic")
    ),

    PGPR = factor(
      PGPR,
      levels = c("No", "Yes")
    ),

    Cultivar = factor(
      Cultivar,
      levels = c("Blueray", "Duke", "Patriot")
    )
  )

# ============================================================
# 2. LONG FORMAT
# ============================================================

long <- X %>%
  pivot_longer(
    cols = -Function,
    names_to = "SampleID",
    values_to = "Relative"
  ) %>%
  left_join(
    meta,
    by = "SampleID"
  )

stopifnot(
  n_distinct(long$SampleID) == 48
)

functions <- X$Function

cat("===== INPUT QC =====\n")
cat("Functions:", length(functions), "\n")
cat("Samples:", n_distinct(long$SampleID), "\n")
cat(
  "Missing metadata:",
  sum(is.na(long$Traitement)),
  "\n"
)

# ============================================================
# 3. DOTPLOT: FERTILISATION x PGPR
# ============================================================

dot <- long %>%
  group_by(
    Function,
    Fertilisation,
    PGPR
  ) %>%
  summarise(
    Mean = mean(Relative),
    SE = sd(Relative) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(
    Group = paste(
      Fertilisation,
      PGPR,
      sep = " / "
    ),
    MeanPercent = 100 * Mean,
    SEPercent = 100 * SE
  )

p1 <- plot_ly(
  dot,
  x = ~MeanPercent,
  y = ~Function,
  color = ~Group,
  type = "scatter",
  mode = "markers",
  error_x = list(
    type = "data",
    array = dot$SEPercent,
    visible = TRUE
  ),
  hovertemplate = paste(
    "<b>%{y}</b><br>",
    "Group: %{text}<br>",
    "Mean: %{x:.4f}%<br>",
    "<extra></extra>"
  ),
  text = ~Group
) %>%
  layout(
    title = list(
      text = "Predicted functional potential by fertilisation and PGPR"
    ),
    xaxis = list(
      title = "Mean relative pathway abundance (%)"
    ),
    yaxis = list(
      title = ""
    ),
    legend = list(
      title = list(
        text = "Treatment"
      )
    )
  )

saveWidget(
  p1,
  file.path(
    outdir,
    "EMP_CoreFunctions_dotplot_Fertilisation_PGPR.html"
  ),
  selfcontained = TRUE
)

# ============================================================
# 4. HEATMAP: TREATMENT x CULTIVAR
# ============================================================

hm <- long %>%
  group_by(
    Function,
    Traitement,
    Cultivar
  ) %>%
  summarise(
    Mean = mean(Relative),
    .groups = "drop"
  ) %>%
  mutate(
    MeanPercent = 100 * Mean,
    Group = paste(
      Traitement,
      Cultivar,
      sep = " | "
    )
  )

group_order <- expand.grid(
  Traitement = c(
    "CTRLMin",
    "PGPRMin",
    "CTROrg",
    "PGPROrg"
  ),
  Cultivar = c(
    "Blueray",
    "Duke",
    "Patriot"
  ),
  KEEP.OUT.ATTRS = FALSE
) %>%
  mutate(
    Group = paste(
      Traitement,
      Cultivar,
      sep = " | "
    )
  ) %>%
  pull(Group)

hm_matrix <- hm %>%
  select(
    Function,
    Group,
    MeanPercent
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = MeanPercent
  )

hm_matrix <- hm_matrix[
  match(functions, hm_matrix$Function),
]

z_raw <- as.matrix(
  hm_matrix[
    ,
    intersect(
      group_order,
      names(hm_matrix)
    ),
    drop = FALSE
  ]
)

rownames(z_raw) <- hm_matrix$Function

p2 <- plot_ly(
  x = colnames(z_raw),
  y = rownames(z_raw),
  z = z_raw,
  type = "heatmap",
  hovertemplate = paste(
    "<b>%{y}</b><br>",
    "%{x}<br>",
    "Relative abundance: %{z:.4f}%<br>",
    "<extra></extra>"
  )
) %>%
  layout(
    title = list(
      text = "Predicted functional potential across treatment × cultivar"
    ),
    xaxis = list(
      title = "",
      tickangle = -45
    ),
    yaxis = list(
      title = ""
    )
  )

saveWidget(
  p2,
  file.path(
    outdir,
    "EMP_CoreFunctions_heatmap_Treat_x_Cult.html"
  ),
  selfcontained = TRUE
)

# ============================================================
# 5. Z-SCORE HEATMAP
# ============================================================

z_scaled <- t(
  scale(
    t(z_raw)
  )
)

p3 <- plot_ly(
  x = colnames(z_scaled),
  y = rownames(z_scaled),
  z = z_scaled,
  type = "heatmap",
  zmid = 0,
  hovertemplate = paste(
    "<b>%{y}</b><br>",
    "%{x}<br>",
    "Z-score: %{z:.3f}<br>",
    "<extra></extra>"
  )
) %>%
  layout(
    title = list(
      text = "Standardized predicted functional potential across treatment × cultivar"
    ),
    xaxis = list(
      title = "",
      tickangle = -45
    ),
    yaxis = list(
      title = ""
    )
  )

saveWidget(
  p3,
  file.path(
    outdir,
    "EMP_CoreFunctions_heatmap_Treat_x_Cult_Zscore.html"
  ),
  selfcontained = TRUE
)

# ============================================================
# 6. FERTILISATION EFFECT PLOT
# ============================================================

fert <- stats %>%
  filter(
    Effect == "Fertilisation"
  ) %>%
  mutate(
    Significant = ifelse(
      FDR < 0.05,
      "FDR < 0.05",
      "NS"
    ),
    Direction = case_when(
      FDR < 0.05 & Estimate > 0 ~ "Higher in Organic",
      FDR < 0.05 & Estimate < 0 ~ "Higher in Mineral",
      TRUE ~ "Not significant"
    ),
    Label = paste0(
      "Estimate = ",
      sprintf("%.3f", Estimate),
      "<br>FDR = ",
      format(
        FDR,
        digits = 3,
        scientific = TRUE
      )
    )
  ) %>%
  arrange(Estimate)

p4 <- plot_ly(
  fert,
  x = ~Estimate,
  y = ~reorder(
    Function,
    Estimate
  ),
  type = "bar",
  orientation = "h",
  color = ~Direction,
  text = ~Label,
  hovertemplate = paste(
    "<b>%{y}</b><br>",
    "%{text}<br>",
    "<extra></extra>"
  )
) %>%
  layout(
    title = list(
      text = "Overall fertilisation effect on targeted predicted functional potential"
    ),
    xaxis = list(
      title = "Model effect: Organic vs Mineral<br>(log₂ stabilized log-ratio scale)",
      zeroline = TRUE
    ),
    yaxis = list(
      title = ""
    ),
    legend = list(
      title = list(
        text = "Direction"
      )
    )
  )

saveWidget(
  p4,
  file.path(
    outdir,
    "EMP_CoreFunctions_Fertilisation_effect.html"
  ),
  selfcontained = TRUE
)

# ============================================================
# 7. SUMMARY
# ============================================================

sigfert <- fert %>%
  filter(FDR < 0.05)

cat("\n===== FINAL EMP HTML QC =====\n")
cat("Functions:", length(functions), "\n")
cat(
  "Significant fertilisation functions:",
  nrow(sigfert),
  "/",
  length(functions),
  "\n"
)

cat(
  "Higher in Organic:",
  sum(
    sigfert$Estimate > 0
  ),
  "\n"
)

cat(
  "Higher in Mineral:",
  sum(
    sigfert$Estimate < 0
  ),
  "\n"
)

cat("\nGenerated files:\n")

for (f in list.files(
  outdir,
  pattern = "\\.html$",
  full.names = TRUE
)) {
  cat(f, "\n")
}
