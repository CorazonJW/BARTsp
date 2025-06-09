#' Construct BART Geneset Input
#'
#' This function prepares input data for BART region analysis by extracting relevant peaks from an expression matrix.
#'
#' @param genes A list containing significant features and their statistics.
#' @return A list containing up-regulated and down-regulated gene sets.
#'
#' @export
construct_BART_geneset_input <- function(genes) {  
  # Input validation
  if (missing(genes)) {
    stop("Error: genes argument is required.")
  }
  
  if (!is.list(genes)) {
    stop("Error: genes must be a list.")
  }
  
  required_elements <- c("significant_features", "statistics", "cor.rho")
  if (!all(required_elements %in% names(genes))) {
    stop("Error: genes must contain 'significant_features', 'statistics', and 'cor.rho'.")
  }
  
  if (length(genes$significant_features) == 0) {
    stop("Error: No significant features found in input.")
  }
  
  if (length(genes$significant_features) != length(genes$statistics) || 
      length(genes$significant_features) != length(genes$cor.rho)) {
    stop("Error: Lengths of significant_features, statistics, and cor.rho must match.")
  }
  
  # Check for missing values
  if (any(is.na(genes$cor.rho))) {
    warning("NA values found in cor.rho. These will be removed.")
    valid_idx <- !is.na(genes$cor.rho)
    genes$significant_features <- genes$significant_features[valid_idx]
    genes$statistics <- genes$statistics[valid_idx]
    genes$cor.rho <- genes$cor.rho[valid_idx]
  }
  
  # increasingly-expressed genes
  up_idx <- genes$cor.rho > 0
  if (sum(up_idx) == 0) {
    warning("No up-regulated genes found.")
  }
  up_genes <- list(
    significant_features = genes$significant_features[up_idx], 
    statistics = genes$statistics[up_idx], 
    cor.rho = genes$cor.rho[up_idx]
  )

  # decreasingly-expressed genes
  down_idx <- genes$cor.rho < 0
  if (sum(down_idx) == 0) {
    warning("No down-regulated genes found.")
  }
  down_genes <- list(
    significant_features = genes$significant_features[down_idx], 
    statistics = genes$statistics[down_idx], 
    cor.rho = genes$cor.rho[down_idx]
  )

  output <- list(up_gene = up_genes, down_gene = down_genes)
  
  return(output)
}

#' Construct BART Region Input
#'
#' This function prepares input data for BART region analysis by extracting relevant peaks from an expression matrix.
#'
#' @param regions A list containing significant regions and their statistics.
#' @param description_up A character string describing up-regulated regions.
#' @param description_down A character string describing down-regulated regions.
#' @return A list containing up-regulated and down-regulated region sets.
#'
#' @export
construct_BART_region_input <- function(regions, description_up, description_down) {  
  # Input validation
  if (missing(regions) || missing(description_up) || missing(description_down)) {
    stop("Error: All arguments are required.")
  }
  
  if (!is.list(regions)) {
    stop("Error: regions must be a list.")
  }
  
  if (!all(c("up_regions", "down_regions") %in% names(regions))) {
    stop("Error: regions must contain 'up_regions' and 'down_regions'.")
  }
  
  if (!is.character(description_up) || !is.character(description_down)) {
    stop("Error: description_up and description_down must be character strings.")
  }
  
  # Validate up_regions
  if (!all(c("significant_features", "statistics") %in% names(regions$up_regions))) {
    stop("Error: up_regions must contain 'significant_features' and 'statistics'.")
  }
  
  # Validate down_regions
  if (!all(c("significant_features", "statistics") %in% names(regions$down_regions))) {
    stop("Error: down_regions must contain 'significant_features' and 'statistics'.")
  }
  
  up_peak <- regions$up_regions$significant_features
  down_peak <- regions$down_regions$significant_features
  
  if (length(up_peak) == 0 && length(down_peak) == 0) {
    stop("Error: No significant regions found in input.")
  }
  
  # Process up-regulated regions
  if (length(up_peak) > 0) {
    peak_info_up <- tryCatch({
      do.call(rbind, strsplit(up_peak, "-"))
    }, error = function(e) {
      stop("Error parsing up-regulated regions: ", e$message)
    })
    
    if (ncol(peak_info_up) != 3) {
      stop("Error: Up-regulated regions must be in format 'chromosome-start-end'.")
    }
    
    value_up <- regions$up_regions$statistics
    if (length(value_up) != length(up_peak)) {
      stop("Error: Length of up-regulated statistics must match number of regions.")
    }
    
    output_df_up <- tryCatch({
      data.frame(
        V1 = peak_info_up[,1],  # Chromosome
        V2 = as.integer(peak_info_up[,2]),  # Start position
        V3 = as.integer(peak_info_up[,3]),  # End position
        V4 = description_up,  # Description
        V5 = value_up,  # Peak statistics
        V6 = "."
      )
    }, error = function(e) {
      stop("Error creating up-regulated data frame: ", e$message)
    })
    
    rownames(output_df_up) <- paste0(output_df_up$V1, "-", output_df_up$V2, "-", output_df_up$V3)
    output_df_up <- output_df_up[match(up_peak, rownames(output_df_up)), , drop=FALSE]
  } else {
    output_df_up <- data.frame()
  }

  # Process down-regulated regions
  if (length(down_peak) > 0) {
    peak_info_down <- tryCatch({
      do.call(rbind, strsplit(down_peak, "-"))
    }, error = function(e) {
      stop("Error parsing down-regulated regions: ", e$message)
    })
    
    if (ncol(peak_info_down) != 3) {
      stop("Error: Down-regulated regions must be in format 'chromosome-start-end'.")
    }
    
    value_down <- regions$down_regions$statistics
    if (length(value_down) != length(down_peak)) {
      stop("Error: Length of down-regulated statistics must match number of regions.")
    }
    
    output_df_down <- tryCatch({
      data.frame(
        V1 = peak_info_down[,1],  # Chromosome
        V2 = as.integer(peak_info_down[,2]),  # Start position
        V3 = as.integer(peak_info_down[,3]),  # End position
        V4 = description_down,  # Description
        V5 = value_down,  # Peak statistics
        V6 = "."
      )
    }, error = function(e) {
      stop("Error creating down-regulated data frame: ", e$message)
    })
    
    rownames(output_df_down) <- paste0(output_df_down$V1, "-", output_df_down$V2, "-", output_df_down$V3)
    output_df_down <- output_df_down[match(down_peak, rownames(output_df_down)), , drop=FALSE]
  } else {
    output_df_down <- data.frame()
  }

  output <- list(up_region = output_df_up, down_region = output_df_down)

  return(output)
}

#' Initialize a BART Project
#'
#' This function initializes a BART project for either gene expression or genomic region data.
#'
#' @param name A character string specifying the analysis name.
#' @param genome A character string specifying the genome version (e.g., "hg38", "mm10").
#' @param data A matrix or data frame containing gene expression or genomic region data.
#' @param type A character string specifying the analysis type ("gene_set" or "region").
#'
#' @return A BART project object.
#'
#' @importFrom BARTsc bart load_bart2
#' @export
bart <- function(name, genome, data, type = "geneset") {
  # Input validation
  if (missing(name) || missing(genome) || missing(data)) {
    stop("Error: name, genome, and data arguments are required.")
  }
  
  if (!is.character(name) || length(name) != 1) {
    stop("Error: name must be a single character string.")
  }
  
  if (!is.character(genome) || length(genome) != 1) {
    stop("Error: genome must be a single character string.")
  }
  
  if (!type %in% c("geneset", "region")) {
    stop("Error: type must be either 'geneset' or 'region'.")
  }
  
  # Load BART data with error handling
  tryCatch({
    BARTsc::load_bart2()
  }, error = function(e) {
    stop("Error loading BART data: ", e$message)
  })
  
  # Initialize BART project with error handling
  tryCatch({
    if (type == "geneset") {
      bart_proj <- BARTsc::bart(name = name, genome = genome, gene_data = data)
    } else {
      bart_proj <- BARTsc::bart(name = name, genome = genome, region_data = data)
    }
  }, error = function(e) {
    stop("Error initializing BART project: ", e$message)
  })
  
  return(bart_proj)
}

#' Run BART Analysis on an Existing Project
#'
#' This function executes the appropriate BART analysis step for a given project.
#'
#' @param bart_proj A BART project object created by \code{initialize_BART_project}.
#' @param type A character string specifying the analysis type ("geneset" or "region").
#'
#' @return A processed BART project object.
#'
#' @importFrom BARTsc run_bart_gene_set run_bart_region
#' @export
run_BART <- function(bart_proj, type = "geneset") {
  # Input validation
  if (missing(bart_proj)) {
    stop("Error: bart_proj argument is required.")
  }
  
  if (!type %in% c("geneset", "region")) {
    stop("Error: type must be either 'geneset' or 'region'.")
  }
  
  # Run BART analysis with error handling
  tryCatch({
    if (type == "geneset") {
      bart_proj <- BARTsc::run_bart_gene_set(bart_proj)
    } else {
      bart_proj <- BARTsc::run_bart_region(bart_proj)
    }
  }, error = function(e) {
    stop("Error running BART analysis: ", e$message)
  })
  
  return(bart_proj)
}

#' Retrieve BART Analysis Results
#'
#' This function extracts the results from a completed BART project.
#'
#' @param bart_proj A BART project object that has been processed by \code{run_BART_step}.
#' @param type A character string specifying the analysis type ("geneset" or "region").
#'
#' @return A data frame containing BART analysis results.
#'
#' @importFrom BARTsc get_bart_result
#' @export
get_BART_results <- function(bart_proj, type = "geneset") {
  # Input validation
  if (missing(bart_proj)) {
    stop("Error: bart_proj argument is required.")
  }
  
  if (!type %in% c("geneset", "region")) {
    stop("Error: type must be either 'geneset' or 'region'.")
  }
  
  # Get BART results with error handling
  tryCatch({
    if (type == "geneset") {
      results <- BARTsc::get_bart_result(bart_proj, "geneset")
    } else {
      results <- BARTsc::get_bart_result(bart_proj, "region")
    }
  }, error = function(e) {
    stop("Error retrieving BART results: ", e$message)
  })
  
  if (is.null(results) || nrow(results) == 0) {
    warning("No results found in BART project.")
  }
  
  return(results)
}

#' Plot BART Analysis Results
#'
#' This function plots the results from a completed BART project.
#'
#' @param BART_results A data frame from \code{get_BART_results}, containing at least columns "TF" and "rank_avg_z_p_a_irwinhall_pvalue".
#' @param TF_of_interest A character vector of transcription factors to highlight.
#' @param cutoff A numeric cutoff value for the Irwin-Hall p-value. Default is 0.1.
#' @param top_n Integer. Number of top-ranked TFs to highlight. Default is 6.
#'
#' @return A ggplot object showing the BART analysis results.
#'
#' @importFrom ggplot2 ggplot geom_point geom_hline labs theme_bw aes
#' @importFrom ggrepel geom_label_repel
#' @export
plot_BART_results <- function(BART_results, TF_of_interest = NULL, cutoff = 0.1, top_n = 6) {
  # Input validation
  if (missing(BART_results)) {
    stop("Error: BART_results argument is required.")
  }
  
  if (!is.data.frame(BART_results)) {
    stop("Error: BART_results must be a data frame.")
  }
  
  required_cols <- c("TF", "rank_avg_z_p_a_irwinhall_pvalue")
  if (!all(required_cols %in% colnames(BART_results))) {
    stop("Error: BART_results must contain columns 'TF' and 'rank_avg_z_p_a_irwinhall_pvalue'")
  }
  
  if (!is.null(TF_of_interest) && !is.character(TF_of_interest)) {
    stop("Error: TF_of_interest must be a character vector or NULL.")
  }
  
  if (!is.numeric(cutoff) || cutoff <= 0 || cutoff > 1) {
    stop("Error: cutoff must be a numeric value between 0 and 1.")
  }
  
  if (!is.numeric(top_n) || top_n < 1) {
    stop("Error: top_n must be a positive integer.")
  }
  
  # Check for missing values
  if (any(is.na(BART_results$rank_avg_z_p_a_irwinhall_pvalue))) {
    warning("NA values found in p-values. These will be removed.")
    BART_results <- BART_results[!is.na(BART_results$rank_avg_z_p_a_irwinhall_pvalue), ]
  }
  
  if (nrow(BART_results) == 0) {
    stop("Error: No valid results to plot.")
  }

  dt <- BART_results %>% dplyr::mutate(rank = seq_len(nrow(.)))

  if (!is.null(TF_of_interest)) {
    missing_TFs <- setdiff(TF_of_interest, dt$TF)
    if (length(missing_TFs) > 0) {
      warning("The following TF(s) were not found in BART database and will not be plotted: ",
              paste(missing_TFs, collapse = ", "))
    }
  }

  top_TF_region <- head(dt$TF, top_n)
  label_TFs <- unique(c(top_TF_region, TF_of_interest))
  y_intercept <- -log10(as.numeric(cutoff))

  # Create plot with error handling
  tryCatch({
    p <- ggplot(dt, aes(x = rank, y = -log10(rank_avg_z_p_a_irwinhall_pvalue))) +
      geom_point(color = "black", size = 1.5, alpha = 0.7) +
      geom_point(data = dt %>% dplyr::filter(TF %in% top_TF_region), 
                 aes(x = rank, y = -log10(rank_avg_z_p_a_irwinhall_pvalue)), 
                 color = "orange", size = 1.5) +
      geom_point(data = dt %>% dplyr::filter(TF %in% TF_of_interest), 
                 aes(x = rank, y = -log10(rank_avg_z_p_a_irwinhall_pvalue)), 
                 color = "red", size = 1.5) +
      ggrepel::geom_label_repel(data = dt %>% dplyr::filter(TF %in% label_TFs), 
                               aes(label = TF), size = 2.5, box.padding = 0.05, 
                               max.overlaps = 25) +
      geom_hline(yintercept = y_intercept, linetype = "dashed", linewidth = 0.5) +
      labs(x = "TF Rank", y = "-log10(p-value)", title = "") +
      theme_bw() + 
      theme(panel.grid.major.x = element_blank(), 
            panel.grid.minor.x = element_blank(), 
            panel.grid.major.y = element_blank(), 
            panel.grid.minor.y = element_blank())
  }, error = function(e) {
    stop("Error creating plot: ", e$message)
  })

  return(p)
} 