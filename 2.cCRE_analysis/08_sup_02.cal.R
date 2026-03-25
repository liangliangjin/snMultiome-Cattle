library(data.table)
read_bed_file <- function(bed_file) {
  bed_df <- fread(bed_file, header = FALSE)
  colnames(bed_df) <- c("chrom", "start", "end", "index")
  return(bed_df)
}
read_score_file <- function(score_file) {
  score_df <- fread(score_file, header = FALSE)
  colnames(score_df) <- c("chrom", "position", "position_plus_one", "score")
  score_df <- score_df[, .(chrom, position, score)] 
  return(score_df)
}
calculate_mean_scores <- function(bed_df, score_df) {
  mean_scores <- numeric(nrow(bed_df))
  for (i in 1:nrow(bed_df)) {
    chrom <- bed_df$chrom[i]
    start <- bed_df$start[i]
    end <- bed_df$end[i]
    region_scores <- score_df[score_df$chrom == chrom & 
                              score_df$position >= start & 
                              score_df$position <= end, ]
    if (nrow(region_scores) > 0) {
      mean_scores[i] <- mean(region_scores$score, na.rm = TRUE)
    } else {
      mean_scores[i] <- NA  # If there are no points in the region
    }
    system(paste("echo 'Processing row", i,'/',nrow(bed_df), "' >>", log_file), intern = TRUE)
  }
  return(mean_scores)
}
save_results <- function(results, bed_df, output_file) {
  result_df <- data.table(index = bed_df$index, mean_score = results)
  fwrite(result_df, file = output_file, sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
bed_file <- paste0(chr, '_region.bed')
score_file <- paste0(chr, '_filter.bed')
output_file <- paste0(chr, '_scores.txt')
log_file <- paste0(chr, '_processing_log.txt')
file.create(log_file)
bed_df <- read_bed_file(bed_file)
score_df <- read_score_file(score_file)
mean_scores <- calculate_mean_scores(bed_df, score_df)
save_results(mean_scores, bed_df, output_file)
cat("Mean results have been saved to", output_file, "\n")