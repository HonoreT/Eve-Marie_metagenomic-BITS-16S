library(readr)
library(dplyr)
library(limma)
library(plotly)
library(htmlwidgets)
library(htmltools)

project <- "/project/166726054/Marie-Eve_Analysis"

matrix_file <- file.path(
  project,
  "08_functional_analysis/BITS/FUNGuild/matrices",
  "all_assigned_guild_relative_total_fungal.tsv"
)

metadata_file <- file.path(
  project,
  "06_exports/BITS_final/Metadata.tsv"
)

outdir <- file.path(
  project,
  "10_figures/BITS_FUNGuild_HTML"
)

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 1. INPUTS
# ============================================================

guild <- read_tsv(
  matrix_file,
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
    ),

    Bloc = factor(Bloc)
  )

# ============================================================
# 2. SAMPLE CONCORDANCE
# ============================================================

sample_cols <- intersect(
  names(guild),
  meta$SampleID
)

meta <- meta %>%
  filter(SampleID %in% sample_cols)

# Keep exact metadata order
sample_cols <- meta$SampleID

X <- as.matrix(
  guild[, sample_cols]
)

rownames(X) <- guild[[1]]

cat("===== INPUT =====\n")
cat("Guilds before prevalence filter:", nrow(X), "\n")
cat("Samples:", ncol(X), "\n")

# ============================================================
# 3. PREVALENCE FILTER
# >= 12 / 48 samples = 25%
# ============================================================

prev <- rowSums(X > 0)

keep <- prev >= 12

Xf <- X[
  keep,
  ,
  drop = FALSE
]

cat(
  "Guilds retained >=12/48 samples:",
  nrow(Xf),
  "\n"
)

# ============================================================
# 4. STABILIZED LOGIT-LIKE TRANSFORMATION
#
# Data are proportions relative to total fungal reads.
# Pseudocount is data-derived:
# half of smallest positive retained abundance.
# ============================================================

minpos <- min(
  Xf[Xf > 0],
  na.rm = TRUE
)

eps <- minpos / 2

Y <- log2(
  (Xf + eps) /
  (1 - Xf + eps)
)

cat(
  "Smallest positive relative abundance:",
  format(minpos, scientific = TRUE),
  "\n"
)

cat(
  "Pseudocount used:",
  format(eps, scientific = TRUE),
  "\n"
)

# ============================================================
# 5. MODEL
# ============================================================

design <- model.matrix(
  ~ Fertilisation + PGPR + Cultivar + Bloc,
  data = meta
)

cat("\n===== DESIGN MATRIX =====\n")
print(colnames(design))

fit <- lmFit(
  Y,
  design
)

fit <- eBayes(fit)

# Explicit coefficients
coef_fert <- "FertilisationOrganic"
coef_pgpr <- "PGPRYes"

if (!coef_fert %in% colnames(design)) {
  stop(
    paste(
      "Coefficient not found:",
      coef_fert
    )
  )
}

if (!coef_pgpr %in% colnames(design)) {
  stop(
    paste(
      "Coefficient not found:",
      coef_pgpr
    )
  )
}

# ============================================================
# 6. RESULTS FUNCTION
# ============================================================

extract_results <- function(
  fit,
  coef_name,
  contrast_name
) {

  tt <- topTable(
    fit,
    coef = coef_name,
    number = Inf,
    sort.by = "none"
  )

  tt$Guild <- rownames(tt)

  tt$Prevalence <- prev[
    tt$Guild
  ]

  result <- tt %>%
    transmute(
      Guild = Guild,
      Effect = logFC,
      PValue = P.Value,
      FDR = adj.P.Val,
      Prevalence = Prevalence,
      NegLog10FDR = -log10(
        pmax(
          FDR,
          1e-300
        )
      ),
      Significant = FDR < 0.05,
      Direction = case_when(
        FDR < 0.05 & Effect > 0 ~ "Positive",
        FDR < 0.05 & Effect < 0 ~ "Negative",
        TRUE ~ "Not significant"
      )
    )

  result
}

fert <- extract_results(
  fit,
  coef_fert,
  "Organic_vs_Mineral"
)

pgpr <- extract_results(
  fit,
  coef_pgpr,
  "PGPR_vs_NoPGPR"
)

# ============================================================
# 7. WRITE FINAL TABLES
# ============================================================

fert_tsv <- file.path(
  outdir,
  "FUNGuild_Guild_Fertilisation_limma_FINAL.tsv"
)

pgpr_tsv <- file.path(
  outdir,
  "FUNGuild_Guild_PGPR_limma_FINAL.tsv"
)

write_tsv(
  fert,
  fert_tsv
)

write_tsv(
  pgpr,
  pgpr_tsv
)

# ============================================================
# 8. INTERACTIVE VOLCANO FUNCTION
# ============================================================

make_volcano <- function(
  df,
  title,
  subtitle,
  positive_text,
  negative_text,
  outfile
) {

  n_total <- nrow(df)

  n_sig <- sum(
    df$FDR < 0.05,
    na.rm = TRUE
  )

  n_pos <- sum(
    df$FDR < 0.05 &
    df$Effect > 0,
    na.rm = TRUE
  )

  n_neg <- sum(
    df$FDR < 0.05 &
    df$Effect < 0,
    na.rm = TRUE
  )

  # ----------------------------------------------------------
  # Labels:
  # up to 12 most statistically informative significant guilds
  # ----------------------------------------------------------

  labels <- df %>%
    filter(FDR < 0.05) %>%
    arrange(
      FDR,
      desc(abs(Effect))
    ) %>%
    slice_head(n = 12)

  # All points
  hover_text <- paste0(
    "<b>", df$Guild, "</b>",
    "<br>Estimated effect: ",
    sprintf("%.3f", df$Effect),
    "<br>P-value: ",
    format(
      df$PValue,
      digits = 3,
      scientific = TRUE
    ),
    "<br>FDR: ",
    format(
      df$FDR,
      digits = 3,
      scientific = TRUE
    ),
    "<br>Prevalence: ",
    df$Prevalence,
    "/48"
  )

  # Base plot
  p <- plot_ly(
    data = df,
    x = ~Effect,
    y = ~NegLog10FDR,
    type = "scatter",
    mode = "markers",
    color = ~Direction,

    colors = c(
      "Negative" = "#377EB8",
      "Not significant" = "#BDBDBD",
      "Positive" = "#D73027"
    ),

    marker = list(
      size = 9,
      opacity = 0.82,
      line = list(
        width = 0.6,
        color = "#444444"
      )
    ),

    text = hover_text,
    hoverinfo = "text"
  )

  # Top significant labels
  if (nrow(labels) > 0) {

    p <- add_trace(
      p,
      data = labels,

      x = ~Effect,
      y = ~NegLog10FDR,

      type = "scatter",
      mode = "markers+text",

      marker = list(
        size = 10,
        symbol = "circle-open",
        color = "black",
        line = list(
          width = 1.4
        )
      ),

      text = ~Guild,

      textposition = "top center",

      textfont = list(
        size = 10,
        color = "black"
      ),

      hoverinfo = "skip",

      showlegend = FALSE,

      inherit = FALSE
    )
  }

  max_y <- max(
    df$NegLog10FDR[
      is.finite(df$NegLog10FDR)
    ],
    na.rm = TRUE
  )

  p <- layout(
    p,

    title = list(
      text = paste0(
        "<b>",
        title,
        "</b><br>",
        "<sup>",
        subtitle,
        "</sup>"
      ),
      x = 0.04
    ),

    xaxis = list(
      title = "Estimated differential effect",
      zeroline = FALSE
    ),

    yaxis = list(
      title = "-log10(FDR)",
      zeroline = FALSE
    ),

    shapes = list(

      # FDR = 0.05
      list(
        type = "line",
        x0 = min(df$Effect, na.rm = TRUE),
        x1 = max(df$Effect, na.rm = TRUE),
        y0 = -log10(0.05),
        y1 = -log10(0.05),
        line = list(
          dash = "dash",
          width = 1,
          color = "#555555"
        )
      ),

      # zero effect
      list(
        type = "line",
        x0 = 0,
        x1 = 0,
        y0 = 0,
        y1 = max_y * 1.05,
        line = list(
          dash = "dot",
          width = 1,
          color = "#555555"
        )
      )
    ),

    legend = list(
      title = list(
        text = "FDR < 0.05"
      ),
      orientation = "h",
      x = 0.02,
      y = -0.18
    ),

    margin = list(
      l = 90,
      r = 70,
      b = 130,
      t = 110
    ),

    hovermode = "closest"
  )

  # ----------------------------------------------------------
  # HTML summary box
  # ----------------------------------------------------------

  summary_box <- tags$div(
    style = paste0(
      "font-family:Arial,Helvetica,sans-serif;",
      "background:#f7f7f7;",
      "border:1px solid #dddddd;",
      "padding:14px 18px;",
      "margin:12px 10px 18px 10px;",
      "line-height:1.55;",
      "font-size:14px;"
    ),

    tags$b(
      "Guild-level differential representation"
    ),

    tags$br(),
    tags$br(),

    paste0(
      "Guilds tested after prevalence filtering: ",
      n_total
    ),

    tags$br(),

    paste0(
      "FDR < 0.05: ",
      n_sig,
      " / ",
      n_total
    ),

    tags$br(),

    paste0(
      positive_text,
      ": ",
      n_pos
    ),

    tags$br(),

    paste0(
      negative_text,
      ": ",
      n_neg
    ),

    tags$br(),
    tags$br(),

    tags$span(
      style = "font-size:13px;color:#555555;",
      paste0(
        "Model: Fertilisation + PGPR + Cultivar + Bloc. ",
        "Values are taxonomically inferred FUNGuild guild ",
        "relative abundances expressed against total fungal reads. ",
        "The x-axis is an estimated model effect after ",
        "stabilized logit-like transformation; it is not a raw-count log2 fold-change."
      )
    )
  )

  p <- prependContent(
    p,
    summary_box
  )

  saveWidget(
    p,
    outfile,
    selfcontained = TRUE
  )

  list(
    total = n_total,
    significant = n_sig,
    positive = n_pos,
    negative = n_neg,
    labels = labels
  )
}

# ============================================================
# 9. FERTILISATION VOLCANO
# ============================================================

fert_html <- file.path(
  outdir,
  "FUNGuild_Guild_Fertilisation_volcano_FINAL.html"
)

fert_stats <- make_volcano(
  fert,

  title =
    "FUNGuild guild-level response to fertilisation",

  subtitle =
    "Positive effect = higher under organic fertilisation; negative effect = higher under mineral fertilisation",

  positive_text =
    "Higher under organic fertilisation",

  negative_text =
    "Higher under mineral fertilisation",

  outfile = fert_html
)

# ============================================================
# 10. PGPR VOLCANO
# ============================================================

pgpr_html <- file.path(
  outdir,
  "FUNGuild_Guild_PGPR_volcano_FINAL.html"
)

pgpr_stats <- make_volcano(
  pgpr,

  title =
    "FUNGuild guild-level response to PGPR",

  subtitle =
    "Positive effect = higher with PGPR; negative effect = higher without PGPR",

  positive_text =
    "Higher with PGPR",

  negative_text =
    "Higher without PGPR",

  outfile = pgpr_html
)

# ============================================================
# 11. CONSOLE SUMMARY
# ============================================================

cat("\n")
cat("========================================\n")
cat("FINAL VOLCANO RESULTS\n")
cat("========================================\n")

cat("\nFERTILISATION\n")
cat(
  "FDR < 0.05:",
  fert_stats$significant,
  "/",
  fert_stats$total,
  "\n"
)
cat(
  "Higher organic:",
  fert_stats$positive,
  "\n"
)
cat(
  "Higher mineral:",
  fert_stats$negative,
  "\n"
)

cat("\nPGPR\n")
cat(
  "FDR < 0.05:",
  pgpr_stats$significant,
  "/",
  pgpr_stats$total,
  "\n"
)
cat(
  "Higher with PGPR:",
  pgpr_stats$positive,
  "\n"
)
cat(
  "Higher without PGPR:",
  pgpr_stats$negative,
  "\n"
)

cat("\n===== TOP LABELLED FERTILISATION GUILDS =====\n")

print(
  fert_stats$labels %>%
    select(
      Guild,
      Effect,
      PValue,
      FDR,
      Prevalence
    )
)

cat("\n===== TOP LABELLED PGPR GUILDS =====\n")

print(
  pgpr_stats$labels %>%
    select(
      Guild,
      Effect,
      PValue,
      FDR,
      Prevalence
    )
)

cat("\n===== OUTPUT FILES =====\n")

cat(
  fert_html,
  "\n"
)

cat(
  pgpr_html,
  "\n"
)

cat(
  fert_tsv,
  "\n"
)

cat(
  pgpr_tsv,
  "\n"
)
