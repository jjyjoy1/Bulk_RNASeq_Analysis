library(dplyr)

# Define the counts directory
counts_dir <- "counts/"

# Get all featureCounts files
files <- list.files(counts_dir, pattern = "_featureCounts.txt$", full.names = TRUE)

# Read each file and extract gene counts
count_list <- lapply(files, function(file) {
  df <- read.delim(file, skip=1, row.names=1)  # Skip first row if it's a header
  df <- df[, c("V7")]  # FeatureCounts stores counts in column 7
  colnames(df) <- gsub("_featureCounts.txt", "", basename(file))  # Sample name
  return(df)
})

# Merge all count tables by gene ID
count_matrix <- Reduce(function(x, y) merge(x, y, by="row.names", all=TRUE), count_list)
rownames(count_matrix) <- count_matrix$Row.names
count_matrix <- count_matrix[, -1]

# Create sample metadata
sample_names <- colnames(count_matrix)
condition <- ifelse(grepl("control", sample_names), "control", "test")  # Adjust based on naming pattern
sample_metadata <- data.frame(sample=sample_names, condition=condition, row.names=sample_names)

# Save the output files
write.table(count_matrix, "combined_counts.txt", sep="\t", quote=FALSE, row.names=TRUE)
write.table(sample_metadata, "sample_metadata.txt", sep="\t", quote=FALSE, row.names=TRUE)

