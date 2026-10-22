library(readr)
library(dplyr)
library(tidyr)
library(plotly)
library(htmlwidgets)
library(htmltools)

project <- "/project/166726054/Marie-Eve_Analysis"

matrix_file <- file.path(
  project,
  "08_functional_analysis/BITS/FUNGuild/targeted_functions",
  "FUNGuild_core_functions_relative_total_fungal.tsv"
)

stats_file <- file.path(
  project,
  "08_functional_analysis/BITS/FUNGuild/targeted_functions",
  "FUNGuild_core_functions_limma.tsv"
)

metadata_file <- file.path(
  project,
  "06_exports/BITS_final/Metadata.tsv"
)

outdir <- file.path(
  project,
  "10_figures/BITS_FUNGuild_HTML"
)

dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

# ============================================================
# 1. INPUTS
# ============================================================

x <- read_tsv(matrix_file, show_col_types=FALSE)
stats <- read_tsv(stats_file, show_col_types=FALSE)

meta <- read_tsv(
  metadata_file,
  comment="#",
  show_col_types=FALSE
) %>%
  mutate(
    SampleID = as.character(`sample-id`),

    Fertilisation = case_when(
      Traitement %in% c("CTRLMin","PGPRMin") ~ "Mineral",
      Traitement %in% c("CTROrg","PGPROrg") ~ "Organic"
    ),

    PGPR = case_when(
      Traitement %in% c("CTRLMin","CTROrg") ~ "No PGPR",
      Traitement %in% c("PGPRMin","PGPROrg") ~ "PGPR"
    )
  )

# ============================================================
# 2. CONCORDANCE
# ============================================================

sample_cols <- intersect(names(x), meta$SampleID)

meta <- meta %>%
  filter(SampleID %in% sample_cols)

sample_cols <- meta$SampleID

X <- as.matrix(
  x[, sample_cols]
)

rownames(X) <- x$Function

cat("===== CONCORDANCE =====\n")
cat("Functions:", nrow(X), "\n")
cat("Samples:", ncol(X), "\n")
cat(
  "Exact concordance:",
  identical(colnames(X), meta$SampleID),
  "\n"
)

# ============================================================
# 3. LONG FORMAT
# ============================================================

long <- as.data.frame(X) %>%
  tibble::rownames_to_column("Function") %>%
  pivot_longer(
    cols=-Function,
    names_to="SampleID",
    values_to="RelativeAbundance"
  ) %>%
  left_join(
    meta %>%
      select(
        SampleID,
        Traitement,
        Fertilisation,
        PGPR,
        Cultivar,
        Bloc
      ),
    by="SampleID"
  )

if (any(is.na(long$Traitement))) {
  stop("Metadata join produced missing treatment values.")
}

# Keep everything as CHARACTER for Plotly
long$Function      <- as.character(long$Function)
long$Traitement    <- as.character(long$Traitement)
long$Cultivar      <- as.character(long$Cultivar)
long$Fertilisation <- as.character(long$Fertilisation)
long$PGPR          <- as.character(long$PGPR)

function_order <- c(
  "Any Saprotroph",
  "Endophyte",
  "Plant Pathogen",
  "Fungal Parasite",
  "Ericoid Mycorrhizal",
  "Ectomycorrhizal"
)

# ============================================================
# 4. DOTPLOT: FERTILISATION x PGPR
# ============================================================

dot <- long %>%
  mutate(
    Group = case_when(
      Fertilisation=="Mineral" & PGPR=="No PGPR" ~ "Mineral | No PGPR",
      Fertilisation=="Mineral" & PGPR=="PGPR"    ~ "Mineral | PGPR",
      Fertilisation=="Organic" & PGPR=="No PGPR" ~ "Organic | No PGPR",
      Fertilisation=="Organic" & PGPR=="PGPR"    ~ "Organic | PGPR"
    )
  ) %>%
  group_by(Function, Group) %>%
  summarise(
    Mean = mean(RelativeAbundance),
    SD   = sd(RelativeAbundance),
    N    = n(),
    .groups="drop"
  ) %>%
  group_by(Function) %>%
  mutate(
    Zscore = ifelse(
      sd(Mean) > 0,
      (Mean - mean(Mean)) / sd(Mean),
      0
    )
  ) %>%
  ungroup()

group_order <- c(
  "Mineral | No PGPR",
  "Mineral | PGPR",
  "Organic | No PGPR",
  "Organic | PGPR"
)

write_tsv(
  dot,
  file.path(
    outdir,
    "FUNGuild_CoreFunctions_Fertilisation_PGPR_mean.tsv"
  )
)

# Explicit numeric marker sizes
min_size <- 12
max_size <- 46

m <- dot$Mean

if (max(m) > min(m)) {
  marker_size <- min_size +
    (m - min(m)) /
    (max(m) - min(m)) *
    (max_size - min_size)
} else {
  marker_size <- rep(25, length(m))
}

hover_dot <- paste0(
  "<b>", dot$Function, "</b>",
  "<br>", dot$Group,
  "<br>Mean relative abundance: ",
  sprintf("%.2f", 100*dot$Mean), "%",
  "<br>SD: ",
  sprintf("%.2f", 100*dot$SD),
  " percentage points",
  "<br>N = ", dot$N,
  "<br>Row Z-score: ",
  sprintf("%.2f", dot$Zscore)
)

p_dot <- plot_ly(
  x = dot$Group,
  y = dot$Function,

  type = "scatter",
  mode = "markers",

  text = hover_dot,
  hoverinfo = "text",

  marker = list(
    size = marker_size,
    color = dot$Zscore,
    colorscale = list(
      c(0, "#2166AC"),
      c(0.5, "#F7F7F7"),
      c(1, "#B2182B")
    ),
    cmin = -max(abs(dot$Zscore), na.rm=TRUE),
    cmax =  max(abs(dot$Zscore), na.rm=TRUE),
    colorbar = list(
      title = "Row<br>Z-score"
    ),
    line = list(
      width = 0.5,
      color = "#444444"
    )
  )
) %>%
  layout(
    title = list(
      text = paste0(
        "<b>Targeted fungal ecological functions</b>",
        "<br><sup>Fertilisation x PGPR</sup>"
      )
    ),

    xaxis = list(
      title="",
      type="category",
      categoryorder="array",
      categoryarray=group_order,
      tickangle=-20
    ),

    yaxis = list(
      title="",
      type="category",
      categoryorder="array",
      categoryarray=rev(function_order),
      automargin=TRUE
    ),

    margin=list(
      l=190,
      r=120,
      b=130,
      t=100
    )
  )

saveWidget(
  p_dot,
  file.path(
    outdir,
    "FUNGuild_CoreFunctions_dotplot_Fertilisation_PGPR.html"
  ),
  selfcontained=TRUE
)

cat("Dotplot: OK\n")

# ============================================================
# 5. TREATMENT x CULTIVAR MEANS
# ============================================================

tc <- long %>%
  group_by(
    Function,
    Traitement,
    Cultivar
  ) %>%
  summarise(
    Mean=mean(RelativeAbundance),
    SD=sd(RelativeAbundance),
    N=n(),
    .groups="drop"
  ) %>%
  mutate(
    Group=paste(
      Traitement,
      Cultivar,
      sep="__"
    )
  )

group_tc_order <- c(
  "CTRLMin__Blueray",
  "CTRLMin__Duke",
  "CTRLMin__Patriot",

  "PGPRMin__Blueray",
  "PGPRMin__Duke",
  "PGPRMin__Patriot",

  "CTROrg__Blueray",
  "CTROrg__Duke",
  "CTROrg__Patriot",

  "PGPROrg__Blueray",
  "PGPROrg__Duke",
  "PGPROrg__Patriot"
)

tc_wide <- tc %>%
  select(Function, Group, Mean) %>%
  pivot_wider(
    names_from=Group,
    values_from=Mean
  )

# Force exact function order
tc_wide <- tc_wide[
  match(
    function_order,
    tc_wide$Function
  ),
]

mat <- as.matrix(
  tc_wide[, group_tc_order]
)

storage.mode(mat) <- "numeric"

rownames(mat) <- tc_wide$Function

mat_pct <- 100 * mat

write_tsv(
  as.data.frame(mat) %>%
    tibble::rownames_to_column("Function"),
  file.path(
    outdir,
    "FUNGuild_CoreFunctions_Treatment_x_Cultivar_mean.tsv"
  )
)

# ============================================================
# 6. RAW HEATMAP
# ============================================================

heat_raw <- plot_ly(
  x=colnames(mat_pct),
  y=rownames(mat_pct),
  z=mat_pct,

  type="heatmap",
  colorscale="Viridis",

  hovertemplate=paste(
    "<b>%{y}</b>",
    "<br>Group: %{x}",
    "<br>Mean relative abundance: %{z:.2f}%",
    "<extra></extra>"
  )
) %>%
  layout(
    title=list(
      text=paste0(
        "<b>Targeted fungal ecological functions</b>",
        "<br><sup>Treatment x Cultivar - mean relative abundance</sup>"
      )
    ),

    xaxis=list(
      title="",
      type="category",
      categoryorder="array",
      categoryarray=group_tc_order,
      tickangle=-45
    ),

    yaxis=list(
      title="",
      type="category",
      categoryorder="array",
      categoryarray=rev(function_order),
      automargin=TRUE
    ),

    margin=list(
      l=190,
      r=80,
      b=180,
      t=100
    )
  )

saveWidget(
  heat_raw,
  file.path(
    outdir,
    "FUNGuild_CoreFunctions_heatmap_Treat_x_Cult.html"
  ),
  selfcontained=TRUE
)

cat("Raw heatmap: OK\n")

# ============================================================
# 7. Z-SCORE HEATMAP
# ============================================================

zmat <- t(
  scale(
    t(mat_pct)
  )
)

zmat[!is.finite(zmat)] <- 0

zlim <- max(
  abs(zmat),
  na.rm=TRUE
)

heat_z <- plot_ly(
  x=colnames(zmat),
  y=rownames(zmat),
  z=zmat,

  type="heatmap",

  colorscale=list(
    c(0, "#2166AC"),
    c(0.5, "#F7F7F7"),
    c(1, "#B2182B")
  ),

  zmin=-zlim,
  zmax=zlim,

  hovertemplate=paste(
    "<b>%{y}</b>",
    "<br>Group: %{x}",
    "<br>Row Z-score: %{z:.2f}",
    "<extra></extra>"
  )
) %>%
  layout(
    title=list(
      text=paste0(
        "<b>Targeted fungal ecological functions</b>",
        "<br><sup>Treatment x Cultivar - row Z-score</sup>"
      )
    ),

    xaxis=list(
      title="",
      type="category",
      categoryorder="array",
      categoryarray=group_tc_order,
      tickangle=-45
    ),

    yaxis=list(
      title="",
      type="category",
      categoryorder="array",
      categoryarray=rev(function_order),
      automargin=TRUE
    ),

    margin=list(
      l=190,
      r=80,
      b=180,
      t=100
    )
  )

saveWidget(
  heat_z,
  file.path(
    outdir,
    "FUNGuild_CoreFunctions_heatmap_Treat_x_Cult_Zscore.html"
  ),
  selfcontained=TRUE
)

cat("Z-score heatmap: OK\n")

# ============================================================
# 8. FERTILISATION EFFECT
# ============================================================

fert <- stats %>%
  filter(
    Effect=="Fertilisation"
  ) %>%
  mutate(
    Direction=case_when(
      FDR < 0.05 & Estimate > 0 ~ "Higher under organic",
      FDR < 0.05 & Estimate < 0 ~ "Higher under mineral",
      TRUE ~ "Not significant"
    )
  )

fert <- fert[
  match(
    function_order,
    fert$Function
  ),
]

write_tsv(
  fert,
  file.path(
    outdir,
    "FUNGuild_CoreFunctions_Fertilisation_effects.tsv"
  )
)

# Explicit colors: no Plotly factor-color mapping
point_colors <- ifelse(
  fert$Direction=="Higher under organic",
  "#D73027",
  ifelse(
    fert$Direction=="Higher under mineral",
    "#377EB8",
    "#BDBDBD"
  )
)

hover_eff <- paste0(
  "<b>", fert$Function, "</b>",
  "<br>Estimated differential effect: ",
  sprintf("%.3f", fert$Estimate),
  "<br>P-value: ",
  format(
    fert$PValue,
    digits=3,
    scientific=TRUE
  ),
  "<br>FDR: ",
  format(
    fert$FDR,
    digits=3,
    scientific=TRUE
  ),
  "<br>", fert$Direction
)

p_eff <- plot_ly(
  x=fert$Estimate,
  y=fert$Function,

  type="scatter",
  mode="markers",

  text=hover_eff,
  hoverinfo="text",

  marker=list(
    size=15,
    color=point_colors,
    line=list(
      width=1,
      color="#333333"
    )
  )
) %>%
  layout(
    title=list(
      text=paste0(
        "<b>Fertilisation effect on targeted fungal ecological functions</b>",
        "<br><sup>Positive = higher under organic; ",
        "negative = higher under mineral</sup>"
      )
    ),

    xaxis=list(
      title="Estimated differential effect",
      zeroline=FALSE
    ),

    yaxis=list(
      title="",
      type="category",
      categoryorder="array",
      categoryarray=rev(function_order),
      automargin=TRUE
    ),

    shapes=list(
      list(
        type="line",
        x0=0,
        x1=0,
        y0=-0.5,
        y1=5.5,
        line=list(
          color="#555555",
          width=1,
          dash="dash"
        )
      )
    ),

    annotations=list(
      list(
        x=min(fert$Estimate),
        y=5.7,
        text=paste0(
          "FDR < 0.05: ",
          sum(fert$FDR < 0.05),
          "/6"
        ),
        showarrow=FALSE,
        xanchor="left"
      )
    ),

    margin=list(
      l=190,
      r=80,
      b=100,
      t=120
    )
  )

summary_box <- tags$div(
  style=paste0(
    "font-family:Arial,Helvetica,sans-serif;",
    "background:#f7f7f7;",
    "border:1px solid #dddddd;",
    "padding:14px 18px;",
    "margin:12px 10px 18px 10px;",
    "line-height:1.55;",
    "font-size:14px;"
  ),

  tags$b(
    "Targeted fungal ecological signatures"
  ),

  tags$br(),
  tags$br(),

  paste0(
    "Significant fertilisation effects (FDR < 0.05): ",
    sum(fert$FDR < 0.05),
    " / 6"
  ),

  tags$br(),

  paste0(
    "Higher under organic fertilisation: ",
    sum(
      fert$FDR < 0.05 &
      fert$Estimate > 0
    )
  ),

  tags$br(),

  paste0(
    "Higher under mineral fertilisation: ",
    sum(
      fert$FDR < 0.05 &
      fert$Estimate < 0
    )
  ),

  tags$br(),
  tags$br(),

  tags$span(
    style="font-size:13px;color:#555555;",

    paste0(
      "Model: Fertilisation * PGPR + Cultivar + Bloc. ",
      "Signatures were taxonomically inferred with FUNGuild ",
      "using Probable and Highly Probable assignments. ",
      "Signatures are non-exclusive: one ASV may contribute ",
      "to multiple ecological functions. ",
      "Relative abundances are expressed against total fungal reads."
    )
  )
)

p_eff <- prependContent(
  p_eff,
  summary_box
)

saveWidget(
  p_eff,
  file.path(
    outdir,
    "FUNGuild_CoreFunctions_Fertilisation_effect.html"
  ),
  selfcontained=TRUE
)

cat("Fertilisation effect plot: OK\n")

# ============================================================
# 9. FINAL QC
# ============================================================

cat("\n===== FINAL TARGETED ITS HTML OUTPUTS =====\n")

files <- c(
  "FUNGuild_CoreFunctions_dotplot_Fertilisation_PGPR.html",
  "FUNGuild_CoreFunctions_heatmap_Treat_x_Cult.html",
  "FUNGuild_CoreFunctions_heatmap_Treat_x_Cult_Zscore.html",
  "FUNGuild_CoreFunctions_Fertilisation_effect.html"
)

for (f in files) {

  path <- file.path(
    outdir,
    f
  )

  cat(
    sprintf(
      "%-65s %s\n",
      f,
      ifelse(
        file.exists(path),
        "OK",
        "MISSING"
      )
    )
  )
}

cat("\n===== FERTILISATION SUMMARY =====\n")

cat(
  "Significant FDR < 0.05:",
  sum(fert$FDR < 0.05),
  "/ 6\n"
)

cat(
  "Higher organic:",
  sum(
    fert$FDR < 0.05 &
    fert$Estimate > 0
  ),
  "\n"
)

cat(
  "Higher mineral:",
  sum(
    fert$FDR < 0.05 &
    fert$Estimate < 0
  ),
  "\n"
)

cat("\n===== DONE =====\n")
