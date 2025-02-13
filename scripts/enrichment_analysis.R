library(gprofiler2)

# Load genes for enrichment
top_genes <- read.delim("top_xai_genes.txt", header=TRUE)$Gene

# GO Enrichment
go_results <- gost(query = top_genes, organism = "hsapiens", sources = "GO")
write.table(go_results$result, "go_enrichment_results.txt", sep="\t", quote=FALSE, row.names=FALSE)

# KEGG Pathway Enrichment
kegg_results <- gost(query = top_genes, organism = "hsapiens", sources = "KEGG")
write.table(kegg_results$result, "kegg_enrichment_results.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Disease Ontology (DO) Enrichment
do_results <- gost(query = top_genes, organism = "hsapiens", sources = "DO")
write.table(do_results$result, "do_enrichment_results.txt", sep="\t", quote=FALSE, row.names=FALSE)

