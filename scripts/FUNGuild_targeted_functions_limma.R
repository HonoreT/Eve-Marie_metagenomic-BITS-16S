library(readr)
library(dplyr)
library(limma)

project <- "/project/166726054/Marie-Eve_Analysis"

matrix_file <- file.path(
  project,
  "08_functional_analysis/BITS/FUNGuild/targeted_functions",
  "FUNGuild_core_functions_relative_total_fungal.tsv"
)

metadata_file <- file.path(
  project,
  "06_exports/BITS_final/Metadata.tsv"
)

outdir <- file.path(
  project,
  "08_functional_analysis/BITS/FUNGuild/targeted_functions"
)

# ------------------------------------------------------------
# 1. Inputs
# ------------------------------------------------------------

x <- read_tsv(
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

# ------------------------------------------------------------
# 2. Matrix / metadata concordance
# ------------------------------------------------------------

sample_cols <- intersect(
  names(x),
  meta$SampleID
)

meta <- meta %>%
  filter(SampleID %in% sample_cols)

sample_cols <- meta$SampleID

X <- as.matrix(
  x[, sample_cols]
)

rownames(X) <- x$Function

cat("===== CONCORDANCE =====\n")
cat("Functions:", nrow(X), "\n")
cat("Samples  :", ncol(X), "\n")
cat(
  "Exact sample concordance:",
  identical(colnames(X), meta$SampleID),
  "\n"
)

# ------------------------------------------------------------
# 3. Zero structure
# ------------------------------------------------------------

cat("\n===== ZERO COUNTS =====\n")

for (i in seq_len(nrow(X))) {
  cat(
    sprintf(
      "%-24s %2d / 48\n",
      rownames(X)[i],
      sum(X[i, ] == 0)
    )
  )
}

# ------------------------------------------------------------
# 4. Data-derived pseudocount
# ------------------------------------------------------------

positive <- X[X > 0]

min_positive <- min(
  positive,
  na.rm = TRUE
)

eps <- min_positive / 2

cat("\n===== TRANSFORMATION =====\n")
cat(
  "Smallest positive abundance:",
  format(min_positive, scientific = TRUE),
  "\n"
)

cat(
  "Pseudocount:",
  format(eps, scientific = TRUE),
  "\n"
)

# stabilized logit-like transformation
Y <- log2(
  (X + eps) /
  (1 - X + eps)
)

cat(
  "NA after transformation:",
  sum(is.na(Y)),
  "\n"
)

cat(
  "Inf after transformation:",
  sum(is.infinite(Y)),
  "\n"
)

# ------------------------------------------------------------
# 5. Factorial model
# ------------------------------------------------------------

design <- model.matrix(
  ~ Fertilisation * PGPR + Cultivar + Bloc,
  data = meta
)

cat("\n===== DESIGN MATRIX =====\n")
print(colnames(design))

cat(
  "Design rank:",
  qr(design)$rank,
  "/",
  ncol(design),
  "\n"
)

fit <- lmFit(
  Y,
  design
)

fit <- eBayes(fit)

# ------------------------------------------------------------
# 6. Extract coefficients
# ------------------------------------------------------------

wanted <- c(
  Fertilisation = "FertilisationOrganic",
  PGPR = "PGPRYes",
  Interaction = "FertilisationOrganic:PGPRYes"
)

all_results <- list()

for (effect_name in names(wanted)) {

  coef_name <- wanted[[effect_name]]

  if (!coef_name %in% colnames(design)) {
    stop(
      paste(
        "Missing coefficient:",
        coef_name
      )
    )
  }

  tt <- topTable(
    fit,
    coef = coef_name,
    number = Inf,
    sort.by = "none"
  )

  tt$Function <- rownames(tt)

  res <- tt %>%
    transmute(
      Effect = effect_name,
      Function = Function,
      Estimate = logFC,
      AveExpr = AveExpr,
      t = t,
      PValue = P.Value,
      FDR = adj.P.Val
    )

  all_results[[effect_name]] <- res
}

results <- bind_rows(
  all_results
)

outfile <- file.path(
  outdir,
  "FUNGuild_core_functions_limma.tsv"
)

write_tsv(
  results,
  outfile
)

# ------------------------------------------------------------
# 7. Console report
# ------------------------------------------------------------

cat("\n===== TARGETED FUNCTION RESULTS =====\n")

for (effect_name in names(wanted)) {

  cat(
    "\n---",
    effect_name,
    "---\n"
  )

  tmp <- results %>%
    filter(
      Effect == effect_name
    ) %>%
    arrange(FDR)

  print(
    tmp %>%
      select(
        Function,
        Estimate,
        PValue,
        FDR
      ),
    n = Inf
  )
}

cat(
  "\nOutput:",
  outfile,
  "\n"
)
