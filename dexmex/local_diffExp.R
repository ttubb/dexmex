#!/usr/bin/env Rscript

# Libraries
suppressPackageStartupMessages({
    library(DESeq2, quietly=TRUE)
    library(argparse, quietly=TRUE)
})

# Parse command line arguments
parser <- ArgumentParser(description='DESeq2 Pipeline Command Line Tool')

parser$add_argument('--outdir', required=TRUE, help='Path to the output directory.')
parser$add_argument('--coldata_path', required=TRUE, help='Path to a table of sample information (as described in http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula)')
parser$add_argument('--counts_path', required=TRUE, help='A .tsv file with feature names in the first column and counts for different conditions in subsequent columns. You can use convert_featureCounts.py to create this from featureCounts output.')
parser$add_argument('--feature_to_mag_path', required=TRUE, help='Feature to MAG data path. You can use feature_to_mag.py to create this.')
parser$add_argument('--reference_level', required=TRUE, help='Reference level for condition comparison (see http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula)')

args <- parser$parse_args()


default_normalization <- function(count_matrix, coldata_path, reference_level) {
    #' Normalizes counts with DESeq2 using a simple design formula of ~condition
    #' 
    #' @param count_data_path A path to a count matrix
    #' @param coldata_path A path to a coldata file
    #' @param reference_level The reference level for the condition
    #' 
    #' @return A normalized count matrix
    #' 
    
    # Load column data
    coldata = read.csv(file=coldata_path, sep='\t', row.names=1)
    coldata[,1] = factor(coldata[,1], levels=unique(coldata[,1]))
    coldata[,1] = relevel(coldata[,1], ref=reference_level)

    # Assuming the first column of coldata is 'condition'
    names(coldata)[1] <- 'condition'

    # Create DESeq2 object
    dds_raw = DESeqDataSetFromMatrix(countData = as.matrix(round(count_matrix)),
                                     colData = coldata,
                                     design = ~ condition)

    # Perform size factor estimation
    dds_norm = estimateSizeFactors(dds_raw)
    
    # Extract normalized counts
    normalized_counts = counts(dds_norm, normalized=TRUE)

    return(normalized_counts)
}


default_deseq_without_normalization <- function(normalized_count_matrix, coldata_path, reference_level) {
    #' Runs with DESeq2 using a simple design formula of ~condition
    #' Circumvents the usual DESeq2 normalization.#
    #' 
    #' @param normalized_count_matrix A normalized count matrix
    #' @param coldata_path A path to a coldata file
    #' @param reference_level The reference level for the condition
    #'  
    #' @return A DESeq2 object
    
    # Load column data
    coldata = read.csv(file=coldata_path, sep='\t', row.names=1)
    coldata[,1] = factor(coldata[,1], levels=unique(coldata[,1]))
    coldata[,1] = relevel(coldata[,1], ref=reference_level)

    # Assuming the first column of coldata is 'condition'
    names(coldata)[1] <- 'condition'

    # Create DESeq2 object
    dds_raw = DESeqDataSetFromMatrix(countData = as.matrix(round(normalized_count_matrix)),
                                     colData = coldata,
                                     design = ~ condition)
    dds_raw$condition = relevel(dds_raw$condition, ref=reference_level)

    # Create a matrix with 1s as normalization factors and use them to overwrite
    # the normalization factors in the DESeq2 object
    normalization_factors = matrix(1,
                                   ncol=ncol(normalized_count_matrix),
                                   nrow=nrow(normalized_count_matrix))
    normalizationFactors(dds_raw) = normalization_factors

    # Run DESeq2
    dds = DESeq(dds_raw)
    
    return(dds)
}


default_shrunken_results <- function(dds, reference_level, shrinkage_type='apeglm') {
    #' Computes shrunken results for all pairwise comparisons against a reference level
    #'
    #' @param dds A DESeq2 object
    #' @param reference_level The reference level for the condition
    #' @param shrinkage_type "apeglm", "ashr", or "normal"
    #'
    #' @return A list of DESeq2 results for each comparison
    
    # Check if the shrinkage type is valid
    if(!(shrinkage_type %in% c("apeglm", "ashr", "normal"))) {
        stop("Invalid shrinkage type. Must be one of 'apeglm', 'ashr', or 'normal'.")
    }

    # Extract the contrasts (need to remove "Intercept")
    contrast_names = resultsNames(dds)
    contrast_names = setdiff(contrast_names, "Intercept")
    print(contrast_names)

    # Iterate over comparison levels and compute shrunken log2 fold changes
    results_list <- list()
    for(contrast in contrast_names) {
        res_shrunk <- lfcShrink(dds, coef=contrast, type=shrinkage_type)
        results_list[[contrast]] <- res_shrunk
    }

    return(results_list)
}


write_normalized_counts_to_tsv <- function(normalized_counts,
                                           feature_to_mag,
                                           directory=outdir) {
    #' Writes a normalized count matrix to a TSV file, including MAG information
    #' 
    #' @param normalized_counts A normalized count matrix with features as rows
    #' @param feature_to_mag A data frame mapping features to MAGs
    #' @param directory The directory to write the TSV file to
    #' 
    #' @return The path to the written TSV file
    
    if(!dir.exists(directory)) {
        dir.create(directory, recursive = TRUE)
    }
    
    # Ensure the row names of normalized_counts are used as a feature_id column
    normalized_counts_df <- data.frame(feature_id=rownames(normalized_counts), normalized_counts)

    # Find & add the matching MAG names for each feature_id
    mag_names <- feature_to_mag$mag_id[match(normalized_counts_df$feature_id, feature_to_mag$feature_id)]
    normalized_counts_df$mag_id <- mag_names

    # Write the data frame, now including MAG information, to a TSV file
    normalized_counts_path = file.path(directory, "normalized_counts.tsv")
    write.table(normalized_counts_df, file=normalized_counts_path, sep="\t", quote=FALSE, row.names=FALSE)
    
    return(normalized_counts_path)
}


write_shrunken_results_to_tsv <- function(shrunken_results, directory=outdir) {
  #' Writes results to .tsv files.
  #' 
  #' @param shrunken_results A list of DESeq2 results objects
  #' @param directory The directory to write the files to
  #'  
  #' @return NULL

  if(!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE)
  }
  
  # Iterate over the shrunken results list
  for(contrast_name in names(shrunken_results)) {
    result = shrunken_results[[contrast_name]]
    
    # Convert and fix colnames
    result_df <- as.data.frame(result)
    result_df$feature_id <- rownames(result_df)
    result_df <- result_df[, c("feature_id", names(result))]
    
    # Write
    file_path = file.path(directory, paste0(contrast_name, "_results.tsv"))
    write.table(result_df, file=file_path, sep="\t", quote=FALSE, row.names=FALSE)
  }
}

run_pipeline <- function(outdir, coldata_path, counts_path, feature_to_mag_path, reference_level) {
  # Read in the mapping of features to MAGs
  feature_to_mag = read.csv(file=feature_to_mag_path, sep='\t', row.names=1)
  mag_ids = unique(feature_to_mag$mag_id)
  # Read in the counts data as a matrix
  counts = as.matrix(read.csv(file=counts_path, sep='\t', row.names=1))
  # Ensure the 'mag_id' column in 'feature_to_mag' is treated as a factor
  feature_to_mag$mag_id <- factor(feature_to_mag$mag_id)
  
  # Create a list to store count matrices split by MAG ID
  mag_counts <- list()
  # Split the counts matrix by MAG ID
  for(mag_id in levels(feature_to_mag$mag_id)) {
    cat("Current mag_id:", mag_id, "\n")
    mag_rows <- rownames(feature_to_mag[feature_to_mag$mag_id == mag_id, , drop=FALSE])
    if(length(mag_rows) > 0 && all(mag_rows %in% rownames(counts))) {
      mag_counts[[as.character(mag_id)]] <- counts[mag_rows, ]
    }
  }


  # Normalize count matrices for each MAG ID
  failed_mag_ids <- character()
  normalized_mag_counts <- lapply(mag_counts, function(count_matrix) {
    tryCatch({
      # Attempt normalization
      default_normalization(count_matrix, coldata_path, reference_level)
    }, error = function(e) {
      # Handle errors by reporting and excluding the problematic MAG ID
      cat("Normalization failed for MAG ID:", names(count_matrix), "- Excluding from further analysis.\n")
      failed_mag_ids <<- c(failed_mag_ids, names(count_matrix)) # Append the failed MAG ID
      NULL # Return NULL to exclude this MAG ID
    })
  })

  # Filter out NULL values from failed normalization attempts
  normalized_mag_counts <- Filter(Negate(is.null), normalized_mag_counts)
  
  # Report issues with normalization
  if(length(failed_mag_ids) > 0) {
    writeLines(failed_mag_ids, file.path(outdir, "failed_mags.txt"))
  }

  ## Normalize count matrices for each MAG ID
  #normalized_mag_counts <- lapply(mag_counts, function(count_matrix) {
  #  default_normalization(count_matrix, coldata_path, reference_level)
  #})

  # Combine all normalized count matrices into one large matrix
  combined_matrix <- do.call(rbind, normalized_mag_counts)

  # Write the combined normalized counts to disk
  write_normalized_counts_to_tsv(combined_matrix, feature_to_mag, outdir)

  # Run DESeq2 analysis without normalization on the combined matrix
  dds = default_deseq_without_normalization(combined_matrix, coldata_path, reference_level)

  # Compute shrunken results for differential expression
  shrunken_results = default_shrunken_results(dds, reference_level)

  # Write the shrunken results to disk
  write_shrunken_results_to_tsv(shrunken_results, outdir)
}


# Call run_pipeline with parsed command line arguments
run_pipeline(args$outdir, args$coldata_path, args$counts_path, args$feature_to_mag_path, args$reference_level)
