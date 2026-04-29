library(data.table)
args <- commandArgs(trailingOnly = TRUE)
usage <- paste(
  "Usage:",
  "Rscript format_cattle_gwas_for_magma.R <input_dir> <output_dir> <trait_metadata.tsv>",
  "",
  "Input files must be named like: res_TRAIT_BREED",
  "Example: res_AS_CHA",
  "",
  "Output columns:",
  "SNP CHR BP A1 A2 P N Z",
  sep = "\n"
)
if (length(args) < 3) {
  cat(usage, "\n")
  quit(status = 1)
}

input_dir <- normalizePath(args[1], winslash = "/", mustWork = TRUE)
output_dir <- args[2]
trait_meta_file <- normalizePath(args[3], winslash = "/", mustWork = TRUE)

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}
output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)

trait_meta <- fread(trait_meta_file)

required_meta_cols <- c("TRAIT", "CANONICAL_TRAIT", "N")
missing_meta_cols <- setdiff(required_meta_cols, names(trait_meta))
if (length(missing_meta_cols) > 0) {
  stop(sprintf(
    "Trait metadata file is missing columns: %s",
    paste(missing_meta_cols, collapse = ", ")
  ))
}

trait_meta[, TRAIT := as.character(TRAIT)]
trait_meta[, CANONICAL_TRAIT := as.character(CANONICAL_TRAIT)]
trait_meta[, N := as.integer(N)]

extract_file_parts <- function(filename) {
  base <- basename(filename)
  m <- regexec("^res_(.+)_([A-Za-z]+)$", base)
  hits <- regmatches(base, m)[[1]]
  if (length(hits) != 3) {
    return(list(trait = NA_character_, breed = NA_character_))
  }
  list(
    trait = hits[2],
    breed = toupper(hits[3])
  )
}

process_one_file <- function(path) {
  file_info <- extract_file_parts(path)
  trait_raw <- file_info$trait
  breed <- file_info$breed

  if (is.na(trait_raw) || is.na(breed)) {
    warning(sprintf("Skip %s: filename does not match res_TRAIT_BREED", basename(path)))
    return(NULL)
  }

  meta_row <- trait_meta[TRAIT == trait_raw]
  if (nrow(meta_row) == 0) {
    warning(sprintf("Skip %s: trait %s not found in metadata", basename(path), trait_raw))
    return(NULL)
  }

  dt <- fread(path)

  required_cols <- c("SNP", "A1", "A2", "Freq", "b", "se", "p")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0) {
    warning(sprintf(
      "Skip %s: missing columns %s",
      basename(path),
      paste(missing_cols, collapse = ", ")
    ))
    return(NULL)
  }

  snp_parts <- tstrsplit(dt$SNP, ":", fixed = TRUE, keep = 1:4)

  dt[, CHR := snp_parts[[1]]]
  dt[, BP := suppressWarnings(as.integer(snp_parts[[2]]))]
  dt[, REF_FROM_SNP := snp_parts[[3]]]
  dt[, ALT_FROM_SNP := snp_parts[[4]]]

  dt[, A1 := as.character(A1)]
  dt[, A2 := as.character(A2)]
  dt[, BETA := suppressWarnings(as.numeric(b))]
  dt[, SE := suppressWarnings(as.numeric(se))]
  dt[, P := suppressWarnings(as.numeric(p))]
  dt[, FREQ := suppressWarnings(as.numeric(Freq))]
  dt[, N := meta_row$N[[1]]]
  dt[, Z := BETA / SE]

  out <- dt[
    !is.na(SNP) &
      !is.na(CHR) &
      !is.na(BP) &
      !is.na(A1) &
      !is.na(A2) &
      !is.na(P) &
      !is.na(N) &
      !is.na(Z) &
      is.finite(Z) &
      !is.na(SE) &
      SE != 0 &
      P > 0 &
      P <= 1,
    .(SNP, CHR, BP, A1, A2, P, N, Z)
  ]

  canonical_trait <- meta_row$CANONICAL_TRAIT[[1]]
  output_file <- file.path(
    output_dir,
    sprintf("%s_%s_magma_input.tsv.gz", canonical_trait, breed)
  )

  fwrite(out, output_file, sep = "\t")

  data.table(
    input_file = basename(path),
    trait_raw = trait_raw,
    canonical_trait = canonical_trait,
    breed = breed,
    sample_size = meta_row$N[[1]],
    total_rows = nrow(dt),
    output_rows = nrow(out),
    output_file = basename(output_file)
  )
}

input_files <- list.files(
  input_dir,
  pattern = "^res_",
  full.names = TRUE
)

if (length(input_files) == 0) {
  stop(sprintf("No files matching '^res_' found in %s", input_dir))
}

summary_list <- lapply(input_files, process_one_file)
summary_dt <- rbindlist(summary_list, fill = TRUE)

summary_file <- file.path(output_dir, "magma_input_summary.tsv")
fwrite(summary_dt, summary_file, sep = "\t")

cat(sprintf("Processed %d file(s)\n", nrow(summary_dt)))
cat(sprintf("Summary file: %s\n", summary_file))