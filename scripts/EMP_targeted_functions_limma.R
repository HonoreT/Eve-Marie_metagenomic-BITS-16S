library(readr)
library(dplyr)
library(limma)

project <- "/project/166726054/Marie-Eve_Analysis"

matrix_file <- file.path(
  project,
  "08_functional_analysis/EMP/targeted_modules",
  "EMP_core_functions_relative.tsv"
)

metadata_file <- file.path(
  project,
  "06_exports/EMP_final/Metadata.tsv"
)

outdir <- file.path(
  project,
  "08_functional_analysis/EMP/targeted_modules"
)

# ============================================================
# 1. INPUTS
# ============================================================

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

# ============================================================
# 2. SAMPLE CONCORDANCE
# ============================================================

sample_cols <- intersect(
  names(x),
  meta$SampleID
)

meta <- meta %>%
  filter(SampleID %in% sample_cols)

# force matrix order to metadata order
sample_cols <- meta$SampleID

X <- as.matrix(
  x[, sample_cols]
)

storage.mode(X) <- "numeric"
rownames(X) <- x$Function

cat("===== CONCORDANCE =====\n")
cat("Functions:", nrow(X), "\n")
cat("Samples  :", ncol(X), "\n")

cat(
  "Exact sample concordance:",
  identical(colnames(X), meta$SampleID),
  "\n"
)

# ============================================================
# 3. ZERO STRUCTURE
# ============================================================

cat("\n===== ZERO COUNTS =====\n")

for (i in seq_len(nrow(X))) {
  cat(
    sprintf(
      "%-32s %2d / %d\n",
      rownames(X)[i],
      sum(X[i, ] == 0),
      ncol(X)
    )
  )
}

# ============================================================
# 4. DATA-DERIVED PSEUDOCOUNT + LOG-RATIO TRANSFORMATION
# ============================================================

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

# log2(module / remainder), stabilized for zeros
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

# ============================================================
# 5. FULL MODEL FOR INTERACTION
# ============================================================

design_full <- model.matrix(
  ~ Fertilisation * PGPR + Cultivar + Bloc,
  data = meta
)

cat("\n===== FULL DESIGN MATRIX =====\n")
print(colnames(design_full))

cat(
  "Design rank:",
  qr(design_full)$rank,
  "/",
  ncol(design_full),
  "\n"
)

fit_full <- lmFit(
  Y,
  design_full
)

fit_full <- eBayes(
  fit_full
)

interaction_coef <- "FertilisationOrganic:PGPRYes"

interaction <- topTable(
  fit_full,
  coef = interaction_coef,
  number = Inf,
  sort.by = "none"
)

interaction$Function <- rownames(interaction)

interaction_out <- interaction %>%
  transmute(
    Effect = "Fertilisation_x_PGPR",
    Function = Function,
    Estimate = logFC,
    AveExpr = AveExpr,
    t = t,
    PValue = P.Value,
    FDR = adj.P.Val
  )

# ============================================================
# 6. ADDITIVE MODEL FOR OVERALL MAIN EFFECTS
# ============================================================

design_add <- model.matrix(
  ~ Fertilisation + PGPR + Cultivar + Bloc,
  data = meta
)

cat("\n===== ADDITIVE DESIGN MATRIX =====\n")
print(colnames(design_add))

cat(
  "Design rank:",
  qr(design_add)$rank,
  "/",
  ncol(design_add),
  "\n"
)

fit_add <- lmFit(
  Y,
  design_add
)

fit_add <- eBayes(
  fit_add
)

wanted <- c(
  Fertilisation = "FertilisationOrganic",
  PGPR = "PGPRYes"
)

main_results <- list()

for (effect_name in names(wanted)) {

  coef_name <- wanted[[effect_name]]

  tt <- topTable(
    fit_add,
    coef = coef_name,
    number = Inf,
    sort.by = "none"
  )

  tt$Function <- rownames(tt)

  main_results[[effect_name]] <- tt %>%
    transmute(
      Effect = effect_name,
      Function = Function,
      Estimate = logFC,
      AveExpr = AveExpr,
      t = t,
      PValue = P.Value,
      FDR = adj.P.Val
    )
}

main_out <- bind_rows(
  main_results
)

# ============================================================
# 7. WRITE OUTPUTS
# ============================================================

write_tsv(
  interaction_out,
  file.path(
    outdir,
    "EMP_core_functions_interaction_limma.tsv"
  )
)

write_tsv(
  main_out,
  file.path(
    outdir,
    "EMP_core_functions_main_effects_limma.tsv"
  )
)

all_out <- bind_rows(
  main_out,
  interaction_out
)

write_tsv(
  all_out,
  file.path(
    outdir,
    "EMP_core_functions_limma_ALL.tsv"
  )
)

# ============================================================
# 8. CONSOLE SUMMARY
# ============================================================

cat("\n===== INTERACTION =====\n")

print(
  as.data.frame(
    interaction_out %>%
      select(
        Function,
        Estimate,
        PValue,
        FDR
      ) %>%
      arrange(FDR)
  ),
  row.names = FALSE
)

cat("\n===== FERTILISATION - OVERALL EFFECT =====\n")

print(
  as.data.frame(
    main_out %>%
      filter(
        Effect == "Fertilisation"
      ) %>%
      select(
        Function,
        Estimate,
        PValue,
        FDR
      ) %>%
      arrange(FDR)
  ),
  row.names = FALSE
)

cat("\n===== PGPR - OVERALL EFFECT =====\n")

print(
  as.data.frame(
    main_out %>%
      filter(
        Effect == "PGPR"
      ) %>%
      select(
        Function,
        Estimate,
        PValue,
        FDR
      ) %>%
      arrange(FDR)
  ),
  row.names = FALSE
)

cat("\n===== SIGNIFICANCE COUNTS =====\n")

cat(
  "Interaction FDR < 0.05:",
  sum(interaction_out$FDR < 0.05),
  "/",
  nrow(interaction_out),
  "\n"
)

cat(
  "Fertilisation FDR < 0.05:",
  sum(
    main_out$Effect == "Fertilisation" &
    main_out$FDR < 0.05
  ),
  "/ 6\n"
)

cat(
  "PGPR FDR < 0.05:",
  sum(
    main_out$Effect == "PGPR" &
    main_out$FDR < 0.05
  ),
  "/ 6\n"
)

cat("\n===== OUTPUTS =====\n")

cat(
  file.path(
    outdir,
    "EMP_core_functions_interaction_limma.tsv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "EMP_core_functions_main_effects_limma.tsv"
  ),
  "\n"
)

cat(
  file.path(
    outdir,
    "EMP_core_functions_limma_ALL.tsv"
  ),
  "\n"
)
