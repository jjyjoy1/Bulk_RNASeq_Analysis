#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_dir = './fastq_dir/'  // Set your input directory
params.output_dir = './output/'    // Set output directory
params.trimmer = 'fastp'           // Choose a trimming tool: fastp, Trimmomatic, Cutadapt
params.aligner = 'STAR'  // Options: 'STAR' or 'HISAT2'
params.genome_index_STAR = './genome/star_index/'
params.genome_index_HISAT2 = './genome/hisat2_index/genome'
params.annotation_gtf = './annotation.gtf'  // GTF annotation file
params.read_counter = 'featureCounts'  // Options: 'featureCounts' or 'HTSeq'

params.threads = 8


process FASTQ_QC {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    path "qc_reports/${sample_id}_fastqc.html", emit: qc_html
    path "qc_reports/${sample_id}_fastqc.zip", emit: qc_zip

    script:
    """
    fastqc -o qc_reports --threads 4 $reads
    """
}

process TRIM_READS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("trimmed/${sample_id}_R1.fastq.gz"), path("trimmed/${sample_id}_R2.fastq.gz")

    script:
    """
    if [ "$params.trimmer" == "fastp" ]; then
        fastp -i ${reads[0]} -I ${reads[1]} -o trimmed/${sample_id}_R1.fastq.gz -O trimmed/${sample_id}_R2.fastq.gz -h trimmed/${sample_id}_report.html -w 4
    else
        trimmomatic PE -threads 4 ${reads[0]} ${reads[1]} \
            trimmed/${sample_id}_R1.fastq.gz trimmed/${sample_id}_R1_unpaired.fastq.gz \
            trimmed/${sample_id}_R2.fastq.gz trimmed/${sample_id}_R2_unpaired.fastq.gz \
            SLIDINGWINDOW:4:20 MINLEN:36
    fi
    """
}


process MULTIQC {
    tag "MultiQC"

    input:
    path qc_reports

    output:
    path "multiqc_report.html", emit: report

    script:
    """
    multiqc qc_reports -o multiqc/
    """
}

process ALIGNMENT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("bam/${sample_id}.bam")

    script:
    """
    if [ "$params.aligner" == "STAR" ]; then
        STAR --runThreadN $params.threads \
             --genomeDir $params.genome_index_STAR \
             --readFilesIn ${reads[0]} ${reads[1]} \
             --readFilesCommand zcat \
             --outFileNamePrefix star_out/${sample_id}_ \
             --outSAMtype BAM SortedByCoordinate
        mv star_out/${sample_id}_Aligned.sortedByCoord.out.bam bam/${sample_id}.bam

    elif [ "$params.aligner" == "HISAT2" ]; then
        hisat2 -p $params.threads \
               -x $params.genome_index_HISAT2 \
               -1 ${reads[0]} -2 ${reads[1]} \
               | samtools sort -@ $params.threads -o bam/${sample_id}.bam
    fi
    """
}

process BAM_QC {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    path "bam_qc/${sample_id}_flagstat.txt", emit: flagstat
    path "bam_qc/${sample_id}_alignment_metrics.txt", emit: metrics

    script:
    """
    samtools flagstat -@ $params.threads $bam > bam_qc/${sample_id}_flagstat.txt
    picard CollectAlignmentSummaryMetrics I=$bam O=bam_qc/${sample_id}_alignment_metrics.txt R=reference.fasta
    """
}

process FEATURECOUNTS {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    path "counts/${sample_id}_featureCounts.txt", emit: counts

    script:
    """
    featureCounts -T $params.threads -a $params.annotation_gtf -o counts/${sample_id}_featureCounts.txt -p -B -C -t exon -g gene_id $bam
    """
}

process HTSEQ_COUNT {
    tag "$sample_id"

    input:
    tuple val(sample_id), path(bam)

    output:
    path "counts/${sample_id}_htseq_counts.txt", emit: counts

    script:
    """
    htseq-count -f bam -r pos -s no -t exon -i gene_id $bam $params.annotation_gtf > counts/${sample_id}_htseq_counts.txt
    """
}

process COMBINE_COUNTS {
    tag "Combine Counts"

    input:
    path count_files

    output:
    path "combined_counts.txt", emit: counts
    path "sample_metadata.txt", emit: metadata

    script:
    """
    Rscript scripts/prepare_counts_matrix.R
    """
}

process DESEQ2 {
    tag "DESeq2"

    input:
    path "combined_counts.txt"
    path "sample_metadata.txt"

    output:
    path "DESeq2_results.txt"
    path "DESeq2_top30_genes.txt"

    script:
    """
    Rscript scripts/deseq2_analysis.R
    """
}

process EDGE_R {
    tag "edgeR"

    input:
    path "combined_counts.txt"
    path "sample_metadata.txt"

    output:
    path "edgeR_results.txt"
    path "edgeR_top30_genes.txt"

    script:
    """
    Rscript scripts/edgeR_analysis.R
    """
}

process LIMMA_VOOM {
    tag "limma/voom"

    input:
    path "combined_counts.txt"
    path "sample_metadata.txt"

    output:
    path "limma_voom_results.txt"
    path "limma_voom_top30_genes.txt"

    script:
    """
    Rscript scripts/limma_voom_analysis.R
    """
}

process ML_DL_ANALYSIS {
    tag "ML/DL Analysis"

    input:
    path "combined_counts.txt"
    path "sample_metadata.txt"

    output:
    path "rf_model.pkl"
    path "xgb_model.pkl"
    path "nn_model.h5"
    path "rf_best_model.pkl"
    path "xgb_best_model.pkl"
    path "top_xai_genes.txt"
    path "gene_enrichment_results.txt"

    script:
    """
    python scripts/ml_dl_analysis_v2.py
    """
}

process GRN_ANALYSIS {
    tag "Gene Regulatory Network"

    input:
    path "selected_genes.txt"
    path "combined_counts.txt"

    output:
    path "gene_correlation_matrix.txt"
    path "gene_correlation_heatmap.png"
    path "gene_regulatory_network_matrix.txt"
    path "gene_regulatory_network.html"

    script:
    """
    python scripts/gene_regulatory_relationship_analysis.py
    """
}

process GENE_ENRICHMENT {
    tag "GO/KEGG/DO Enrichment"

    input:
    path "top_xai_genes.txt"

    output:
    path "go_enrichment_results.txt"
    path "kegg_enrichment_results.txt"
    path "do_enrichment_results.txt"

    script:
    """
    Rscript scripts/enrichment_analysis.R
    """
}

process WGCNA {
    tag "Weighted Gene Co-expression Network Analysis"

    input:
    path "combined_counts.txt"
    path "sample_metadata.txt"

    output:
    path "module_trait_relationships.txt"
    path "wgcna_edges.txt"
    path "wgcna_nodes.txt"

    script:
    """
    Rscript scripts/weighted_gene_coexpression_network_analysis.R
    """
}


workflow {
    Channel
        .fromPath("$params.input_dir/*.fastq.gz")
        .groupTuple(2, by: { it.split('/')[-1].split('_')[2] }) 
        .map { sample, reads -> tuple(sample, reads.sort()) } 
    
    | FASTQ_QC
    | MULTIQC
    | TRIM_READS
    | ALIGNMENT
    | BAM_QC
    | FEATURECOUNTS

    | COMBINE_COUNTS
    | DESEQ2
    | EDGE_R
    | LIMMA_VOOM
    | ML_DL_ANALYSIS
    | GRN_ANALYSIS
    | GENE_ENRICHMENT
    | WGCNA
}



