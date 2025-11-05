#!/usr/bin/env Rscript
################################################################################
# Travel Ratio Computation Pipeline
#
# Authors: Xiaojun Max Xu, Jinting Li (Original)
#          Refactored and improved by Claude
# Date: November 2025
# Version: 2.0
#
# Description:
#   This script computes the Travel Ratio (TR) from ChIP-seq data.
#   TR is the ratio of gene body signal to TSR (Transcription Start Region)
#   signal, indicating transcriptional elongation efficiency.
#
# Usage:
#   Rscript compute_travel_ratio.R [--config config.yaml]
#
# Requirements:
#   - R packages: yaml, dplyr, tidyr, ggplot2, optparse, logger
#   - External tools: bedtools, bigWigAverageOverBed
################################################################################

#===============================================================================
# SETUP AND LIBRARIES
#===============================================================================

# Required packages
required_packages <- c("yaml", "dplyr", "tidyr", "ggplot2", "optparse",
                       "logger", "readr", "stringr")

# Install missing packages
missing_packages <- required_packages[!required_packages %in%
                                      installed.packages()[,"Package"]]
if(length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# Load libraries
suppressPackageStartupMessages({
  library(yaml)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(optparse)
  library(logger)
  library(readr)
  library(stringr)
})

#===============================================================================
# COMMAND LINE ARGUMENTS
#===============================================================================

option_list <- list(
  make_option(c("-c", "--config"), type = "character", default = "config.yaml",
              help = "Path to configuration YAML file [default: %default]"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output directory (overrides config file)"),
  make_option(c("-v", "--verbose"), action = "store_true", default = FALSE,
              help = "Print detailed progress messages")
)

opt_parser <- OptionParser(option_list = option_list,
                          description = "Compute Travel Ratio from ChIP-seq data")
opt <- parse_args(opt_parser)

#===============================================================================
# CONFIGURATION AND VALIDATION
#===============================================================================

#' Load and validate configuration
#'
#' @param config_file Path to YAML configuration file
#' @return List containing validated configuration parameters
load_config <- function(config_file) {
  log_info("Loading configuration from: {config_file}")

  if (!file.exists(config_file)) {
    log_error("Configuration file not found: {config_file}")
    stop("Configuration file not found")
  }

  config <- yaml::read_yaml(config_file)
  log_success("Configuration loaded successfully")

  return(config)
}

#' Validate tool availability
#'
#' @param tool_name Name or path of the tool
#' @return TRUE if tool is available, FALSE otherwise
validate_tool <- function(tool_name) {
  result <- suppressWarnings(system2("which", tool_name,
                                     stdout = TRUE, stderr = TRUE))

  if (length(result) == 0 || attr(result, "status") == 1) {
    log_warn("Tool not found in PATH: {tool_name}")
    return(FALSE)
  }

  log_debug("Tool found: {tool_name} at {result}")
  return(TRUE)
}

#' Validate all required tools
#'
#' @param config Configuration list
#' @return TRUE if all tools available, stops execution otherwise
validate_tools <- function(config) {
  log_info("Validating external tools...")

  tools_to_check <- c(
    config$tools$bedtools,
    config$tools$bigwig_average_over_bed
  )

  tools_status <- sapply(tools_to_check, validate_tool)

  if (!all(tools_status)) {
    missing_tools <- names(tools_status)[!tools_status]
    log_error("Missing required tools: {paste(missing_tools, collapse = ', ')}")
    stop("Please install missing tools or update paths in config.yaml")
  }

  log_success("All required tools are available")
  return(TRUE)
}

#' Setup output directory structure
#'
#' @param config Configuration list
#' @return Path to output directory
setup_output_directory <- function(config) {
  output_dir <- config$output$output_dir

  if (!dir.exists(output_dir)) {
    log_info("Creating output directory: {output_dir}")
    dir.create(output_dir, recursive = TRUE)
  }

  # Create subdirectories
  subdirs <- c("bed_files", "signals", "plots", "logs")
  for (subdir in subdirs) {
    subdir_path <- file.path(output_dir, subdir)
    if (!dir.exists(subdir_path)) {
      dir.create(subdir_path, recursive = TRUE)
    }
  }

  log_success("Output directory structure created")
  return(output_dir)
}

#' Setup logging configuration
#'
#' @param config Configuration list
setup_logging <- function(config) {
  log_threshold <- switch(
    config$logging$log_level,
    "DEBUG" = DEBUG,
    "INFO" = INFO,
    "WARN" = WARN,
    "ERROR" = ERROR,
    INFO  # default
  )

  log_threshold(log_threshold)

  if (config$logging$save_log) {
    log_file <- file.path(config$output$output_dir, "logs",
                         config$logging$log_file)
    log_appender(appender_tee(log_file))
  }

  log_info("Logging initialized at level: {config$logging$log_level}")
}

#===============================================================================
# CORE ANALYSIS FUNCTIONS
#===============================================================================

#' Convert HOMER peak file to standard BED format and filter by length
#'
#' @param homer_file Path to HOMER peak file
#' @param output_file Path to output BED file
#' @param min_length Minimum gene length in base pairs
#' @return Path to output BED file
convert_and_filter_peaks <- function(homer_file, output_file, min_length = 1000) {
  log_info("Converting HOMER peaks and filtering genes >= {min_length}bp")

  if (!file.exists(homer_file)) {
    log_error("Input file not found: {homer_file}")
    stop("Input file not found")
  }

  # Read HOMER file and convert to BED format
  # HOMER format: name, chr, start, end, strand, ...
  # BED format: chr, start, end, name, score, strand

  cmd <- sprintf(
    "awk -v OFS='\\t' '{print $2, $3, $4, $1, 0, $5}' %s | sort -k1,1 -k2,2n | awk -v OFS='\\t' '$3-$2 > %d {print $0}' > %s",
    homer_file, min_length, output_file
  )

  result <- system(cmd, intern = TRUE)

  if (!file.exists(output_file)) {
    log_error("Failed to create filtered BED file")
    stop("BED file creation failed")
  }

  # Count genes
  n_genes <- as.numeric(system(sprintf("wc -l < %s", output_file), intern = TRUE))
  log_success("Created BED file with {n_genes} genes (length >= {min_length}bp)")

  return(output_file)
}

#' Generate TSR (Transcription Start Region) BED file
#'
#' @param bed_file Path to gene BED file
#' @param output_file Path to output TSR BED file
#' @param upstream Extension upstream of TSS (bp)
#' @param downstream Extension downstream of TSS (bp)
#' @return Path to TSR BED file
generate_tsr_bed <- function(bed_file, output_file, upstream = 50, downstream = 300) {
  log_info("Generating TSR regions (upstream: {upstream}bp, downstream: {downstream}bp)")

  # Read BED file
  bed_data <- read_tsv(bed_file, col_names = c("chr", "start", "end", "name",
                                                "score", "strand"),
                      col_types = "ciicic", show_col_types = FALSE)

  # Calculate TSR coordinates based on strand
  tsr_data <- bed_data %>%
    mutate(
      tsr_start = ifelse(strand == "+", start - upstream, end - downstream),
      tsr_end = ifelse(strand == "+", start + downstream, end + upstream)
    ) %>%
    select(chr, tsr_start, tsr_end, name) %>%
    rename(start = tsr_start, end = tsr_end)

  # Write TSR BED file
  write_tsv(tsr_data, output_file, col_names = FALSE)

  log_success("Created TSR BED file with {nrow(tsr_data)} regions")
  return(output_file)
}

#' Compute signal over BED regions using bigWigAverageOverBed
#'
#' @param bigwig_file Path to bigWig file
#' @param bed_file Path to BED file
#' @param output_file Path to output signal file
#' @param tool_path Path to bigWigAverageOverBed tool
#' @return Path to signal output file
compute_signal_over_bed <- function(bigwig_file, bed_file, output_file,
                                   tool_path = "bigWigAverageOverBed") {
  log_info("Computing signal over BED regions")

  if (!file.exists(bigwig_file)) {
    log_error("BigWig file not found: {bigwig_file}")
    stop("BigWig file not found")
  }

  if (!file.exists(bed_file)) {
    log_error("BED file not found: {bed_file}")
    stop("BED file not found")
  }

  cmd <- sprintf("%s %s %s %s", tool_path, bigwig_file, bed_file, output_file)

  result <- system(cmd, intern = FALSE)

  if (result != 0 || !file.exists(output_file)) {
    log_error("Failed to compute signal over BED regions")
    stop("Signal computation failed")
  }

  log_success("Signal computation completed: {output_file}")
  return(output_file)
}

#' Plot cumulative distribution function of signal values
#'
#' @param signal_file Path to signal file from bigWigAverageOverBed
#' @param output_plot Path to output plot file
#' @param title Plot title
#' @param config Configuration list for plot parameters
plot_signal_cdf <- function(signal_file, output_plot, title = "Signal Distribution",
                           config = NULL) {
  log_info("Generating CDF plot")

  # Read signal data
  signal_data <- read_tsv(signal_file,
                         col_names = c("name", "size", "covered", "sum",
                                      "mean_overall", "mean_exon"),
                         col_types = "ciiddd", show_col_types = FALSE)

  # Create CDF plot
  p <- ggplot(signal_data, aes(x = mean_overall)) +
    stat_ecdf(geom = "point", size = 0.5, color = "#2E86AB") +
    labs(
      title = title,
      x = "Mean Signal",
      y = "Cumulative Proportion",
      caption = sprintf("n = %d regions", nrow(signal_data))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  # Save plot
  if (!is.null(config)) {
    ggsave(output_plot, p,
           width = config$plotting$plot_width,
           height = config$plotting$plot_height,
           dpi = config$plotting$plot_dpi)
  } else {
    ggsave(output_plot, p, width = 8, height = 6, dpi = 300)
  }

  log_success("CDF plot saved: {output_plot}")

  # Return summary statistics
  summary_stats <- summary(signal_data$mean_overall)
  log_info("Signal summary: Min={round(summary_stats[1], 3)}, Median={round(summary_stats[3], 3)}, Mean={round(summary_stats[4], 3)}, Max={round(summary_stats[6], 3)}")

  return(p)
}

#' Select unique transcripts based on signal intensity
#'
#' This function replaces the Python code in the original script.
#' For each gene, it selects the transcript with the highest signal
#' that meets the minimum cutoff threshold.
#'
#' @param signal_file Path to signal file
#' @param output_file Path to output file with unique genes
#' @param cutoff Minimum mean signal threshold
#' @return Path to output file
select_unique_transcripts <- function(signal_file, output_file, cutoff = 0.25) {
  log_info("Selecting unique transcripts with signal >= {cutoff}")

  # Read signal data
  signal_data <- read_tsv(signal_file,
                         col_names = c("region_name", "size", "covered", "sum",
                                      "mean_overall", "mean_exon"),
                         col_types = "ciiddd", show_col_types = FALSE)

  # Extract gene name from region name (format: transcript!gene!...)
  signal_data <- signal_data %>%
    mutate(gene = str_split_fixed(region_name, "!", 3)[, 2])

  # Filter by cutoff and select transcript with highest signal per gene
  unique_transcripts <- signal_data %>%
    filter(mean_overall >= cutoff) %>%
    group_by(gene) %>%
    slice_max(order_by = mean_overall, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Write output in same format as input
  unique_transcripts %>%
    select(region_name, size, covered, sum, mean_overall, mean_exon) %>%
    write_tsv(output_file, col_names = FALSE)

  n_total <- length(unique(signal_data$gene))
  n_selected <- nrow(unique_transcripts)
  n_filtered <- n_total - n_selected

  log_success("Selected {n_selected} unique transcripts from {n_total} genes")
  log_info("Filtered out {n_filtered} genes below cutoff")

  return(output_file)
}

#' Extract BED regions for selected transcripts
#'
#' @param signal_file Path to unique transcripts signal file
#' @param bed_file Path to original BED file
#' @param output_file Path to output BED file
#' @return Path to output BED file
extract_selected_bed_regions <- function(signal_file, bed_file, output_file) {
  log_info("Extracting BED regions for selected transcripts")

  # Read unique transcript names
  unique_data <- read_tsv(signal_file,
                         col_names = c("region_name", "size", "covered",
                                      "sum", "mean_overall", "mean_exon"),
                         col_types = "ciiddd", show_col_types = FALSE)

  # Read BED file
  bed_data <- read_tsv(bed_file,
                      col_names = c("chr", "start", "end", "name", "score", "strand"),
                      col_types = "ciicic", show_col_types = FALSE)

  # Filter BED file for selected transcripts
  selected_bed <- bed_data %>%
    filter(name %in% unique_data$region_name)

  # Remove duplicate genomic coordinates
  unique_bed <- selected_bed %>%
    distinct(chr, start, end, .keep_all = TRUE)

  # Write output
  write_tsv(unique_bed, output_file, col_names = FALSE)

  n_removed <- nrow(selected_bed) - nrow(unique_bed)
  log_success("Extracted {nrow(unique_bed)} unique genomic regions")
  if (n_removed > 0) {
    log_info("Removed {n_removed} duplicate genomic coordinates")
  }

  return(output_file)
}

#' Remove overlapping genes using bedtools
#'
#' @param bed_file Path to input BED file
#' @param output_file Path to output non-overlapping BED file
#' @param distance_cutoff Minimum distance between genes (bp)
#' @param bedtools_path Path to bedtools executable
#' @return Path to non-overlapping BED file
remove_overlapping_genes <- function(bed_file, output_file, distance_cutoff = 1000,
                                    bedtools_path = "bedtools") {
  log_info("Removing overlapping genes (distance cutoff: {distance_cutoff}bp)")

  temp_dir <- dirname(output_file)
  counted_file <- file.path(temp_dir, "temp_counted.bed")
  filtered_file <- file.path(temp_dir, "temp_filtered.bed")

  # Count overlaps
  merge_cmd <- sprintf(
    "%s merge -i %s -d %d -c 1 -o count > %s",
    bedtools_path, bed_file, distance_cutoff, counted_file
  )
  system(merge_cmd)

  # Filter only non-overlapping regions
  filter_cmd <- sprintf("awk '/\\t1$/{print}' %s > %s", counted_file, filtered_file)
  system(filter_cmd)

  # Intersect with original to get full BED records
  intersect_cmd <- sprintf(
    "%s intersect -a %s -b %s -wa > %s",
    bedtools_path, bed_file, filtered_file, output_file
  )
  system(intersect_cmd)

  # Cleanup
  file.remove(counted_file, filtered_file)

  # Count results
  n_input <- as.numeric(system(sprintf("wc -l < %s", bed_file), intern = TRUE))
  n_output <- as.numeric(system(sprintf("wc -l < %s", output_file), intern = TRUE))
  n_removed <- n_input - n_output

  log_success("Removed {n_removed} overlapping genes, retained {n_output} genes")

  return(output_file)
}

#' Generate gene body BED file with extension
#'
#' @param bed_file Path to gene BED file
#' @param tsr_output Path to output TSR BED file
#' @param gb_output Path to output gene body BED file
#' @param tsr_upstream TSR upstream extension (bp)
#' @param tsr_downstream TSR downstream extension (bp)
#' @param gb_extension Gene body extension (bp)
generate_tsr_and_genebody_beds <- function(bed_file, tsr_output, gb_output,
                                          tsr_upstream = 50,
                                          tsr_downstream = 300,
                                          gb_extension = 500) {
  log_info("Generating TSR and gene body BED files")
  log_info("TSR: -{tsr_upstream}bp to +{tsr_downstream}bp from TSS")
  log_info("Gene body: +{tsr_downstream + 1}bp from TSS to TES + {gb_extension}bp")

  # Read BED file
  bed_data <- read_tsv(bed_file,
                      col_names = c("chr", "start", "end", "name", "score", "strand"),
                      col_types = "ciicic", show_col_types = FALSE)

  # Generate TSR and gene body coordinates
  regions <- bed_data %>%
    mutate(
      # TSR coordinates
      tsr_start = ifelse(strand == "+", start - tsr_upstream, end - tsr_downstream),
      tsr_end = ifelse(strand == "+", start + tsr_downstream, end + tsr_upstream),

      # Gene body coordinates
      gb_start = ifelse(strand == "+", start + tsr_downstream + 1, start - gb_extension),
      gb_end = ifelse(strand == "+", end + gb_extension, end - tsr_downstream - 1)
    )

  # Write TSR BED
  tsr_bed <- regions %>%
    select(chr, tsr_start, tsr_end, name) %>%
    rename(start = tsr_start, end = tsr_end)
  write_tsv(tsr_bed, tsr_output, col_names = FALSE)

  # Write gene body BED
  gb_bed <- regions %>%
    select(chr, gb_start, gb_end, name) %>%
    rename(start = gb_start, end = gb_end)
  write_tsv(gb_bed, gb_output, col_names = FALSE)

  log_success("Created TSR and gene body BED files for {nrow(bed_data)} genes")

  return(list(tsr = tsr_output, gb = gb_output))
}

#' Compute Travel Ratio from TSR and gene body signals
#'
#' @param tsr_signal_file Path to TSR signal file
#' @param gb_signal_file Path to gene body signal file
#' @param output_file Path to output TR file
#' @return Data frame with Travel Ratio values
compute_travel_ratio <- function(tsr_signal_file, gb_signal_file, output_file) {
  log_info("Computing Travel Ratio (TR = GB_signal / TSR_signal)")

  # Column names for bigWigAverageOverBed output
  col_names <- c("name", "size", "covered", "sum", "mean_overall", "mean_exon")

  # Read signal files
  tsr_signal <- read_tsv(tsr_signal_file, col_names = col_names,
                        col_types = "ciiddd", show_col_types = FALSE)
  gb_signal <- read_tsv(gb_signal_file, col_names = col_names,
                       col_types = "ciiddd", show_col_types = FALSE)

  # Merge and compute TR
  tr_data <- inner_join(
    tsr_signal %>% select(name, tsr_mean = mean_overall),
    gb_signal %>% select(name, gb_mean = mean_overall),
    by = "name"
  ) %>%
    mutate(
      travel_ratio = gb_mean / tsr_mean,
      # Handle edge cases
      travel_ratio = ifelse(is.infinite(travel_ratio) | is.nan(travel_ratio),
                           NA, travel_ratio)
    )

  # Write output
  tr_output <- tr_data %>%
    select(name, travel_ratio) %>%
    filter(!is.na(travel_ratio))

  write_tsv(tr_output, output_file, col_names = FALSE)

  # Summary statistics
  tr_summary <- summary(tr_output$travel_ratio)
  log_success("Computed TR for {nrow(tr_output)} genes")
  log_info("TR summary: Min={round(tr_summary[1], 3)}, Median={round(tr_summary[3], 3)}, Mean={round(tr_summary[4], 3)}, Max={round(tr_summary[6], 3)}")

  return(tr_data)
}

#' Plot Travel Ratio distribution
#'
#' @param tr_file Path to TR file
#' @param output_plot Path to output plot
#' @param config Configuration list
plot_travel_ratio_cdf <- function(tr_file, output_plot, config = NULL) {
  log_info("Generating Travel Ratio CDF plot")

  # Read TR data
  tr_data <- read_tsv(tr_file, col_names = c("name", "travel_ratio"),
                     col_types = "cd", show_col_types = FALSE)

  # Create CDF plot
  p <- ggplot(tr_data, aes(x = travel_ratio)) +
    stat_ecdf(geom = "point", size = 0.5, color = "#A23B72") +
    labs(
      title = "Travel Ratio Distribution",
      x = "Travel Ratio (Gene Body / TSR)",
      y = "Cumulative Proportion",
      caption = sprintf("n = %d genes", nrow(tr_data))
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )

  # Save plot
  if (!is.null(config)) {
    ggsave(output_plot, p,
           width = config$plotting$plot_width,
           height = config$plotting$plot_height,
           dpi = config$plotting$plot_dpi)
  } else {
    ggsave(output_plot, p, width = 8, height = 6, dpi = 300)
  }

  log_success("Travel Ratio CDF plot saved: {output_plot}")

  return(p)
}

#===============================================================================
# MAIN PIPELINE
#===============================================================================

#' Run complete Travel Ratio computation pipeline
#'
#' @param config_file Path to configuration YAML file
#' @return List containing paths to output files and results
run_travel_ratio_pipeline <- function(config_file = "config.yaml") {

  log_info("=" %R% strrep("=", 70))
  log_info("Starting Travel Ratio Computation Pipeline")
  log_info("=" %R% strrep("=", 70))

  start_time <- Sys.time()

  # Load configuration
  config <- load_config(config_file)

  # Setup output directory
  output_dir <- setup_output_directory(config)

  # Setup logging
  setup_logging(config)

  # Validate tools
  validate_tools(config)

  # Define file paths
  prefix <- config$output$prefix
  bed_dir <- file.path(output_dir, "bed_files")
  signal_dir <- file.path(output_dir, "signals")
  plot_dir <- file.path(output_dir, "plots")

  # Step 1: Convert and filter peaks
  log_info("\n### Step 1: Converting and filtering gene annotations ###")
  bed_filtered <- file.path(bed_dir, paste0(prefix, "_genes_gt",
                           config$gene_filters$min_gene_length, "bp.bed"))
  convert_and_filter_peaks(
    config$input$homer_peak_file,
    bed_filtered,
    config$gene_filters$min_gene_length
  )

  # Step 2: Generate TSR BED for transcript selection
  log_info("\n### Step 2: Generating TSR regions for transcript selection ###")
  tsr_for_selection <- file.path(bed_dir, paste0(prefix, "_tsr_for_selection.bed"))
  generate_tsr_bed(
    bed_filtered,
    tsr_for_selection,
    config$tsr_parameters$upstream_extension,
    config$tsr_parameters$downstream_extension
  )

  # Step 3: Compute signal over genes for selection
  log_info("\n### Step 3: Computing signal for transcript selection ###")
  selection_signal <- file.path(signal_dir, paste0(prefix, "_selection_signal.txt"))
  compute_signal_over_bed(
    config$input$selection_bigwig,
    bed_filtered,  # Use full gene for selection
    selection_signal,
    config$tools$bigwig_average_over_bed
  )

  # Plot signal distribution
  if (config$plotting$generate_plots) {
    selection_plot <- file.path(plot_dir, paste0(prefix, "_selection_signal_cdf.",
                                                 config$plotting$plot_format))
    plot_signal_cdf(selection_signal, selection_plot,
                   "Signal Distribution for Transcript Selection", config)
  }

  # Step 4: Select unique transcripts
  log_info("\n### Step 4: Selecting unique transcripts ###")
  unique_transcripts <- file.path(signal_dir, paste0(prefix, "_unique_transcripts.txt"))
  select_unique_transcripts(
    selection_signal,
    unique_transcripts,
    config$transcript_selection$signal_cutoff
  )

  # Step 5: Extract BED regions for selected transcripts
  log_info("\n### Step 5: Extracting BED regions for selected transcripts ###")
  bed_unique <- file.path(bed_dir, paste0(prefix, "_unique_genes.bed"))
  extract_selected_bed_regions(unique_transcripts, bed_filtered, bed_unique)

  # Step 6: Remove overlapping genes
  log_info("\n### Step 6: Removing overlapping genes ###")
  bed_nonoverlap <- file.path(bed_dir, paste0(prefix, "_nonoverlapping_genes.bed"))
  remove_overlapping_genes(
    bed_unique,
    bed_nonoverlap,
    config$gene_filters$distance_cutoff,
    config$tools$bedtools
  )

  # Step 7: Generate final TSR and gene body BED files
  log_info("\n### Step 7: Generating final TSR and gene body regions ###")
  tsr_final <- file.path(bed_dir, paste0(prefix, "_tsr_final.bed"))
  gb_final <- file.path(bed_dir, paste0(prefix, "_genebody_final.bed"))

  # Calculate gene body extension
  if (config$genebody_parameters$auto_calculate_extension) {
    gb_extension <- config$gene_filters$distance_cutoff / 2
    log_info("Auto-calculated gene body extension: {gb_extension}bp")
  } else {
    gb_extension <- config$genebody_parameters$manual_extension
    log_info("Using manual gene body extension: {gb_extension}bp")
  }

  generate_tsr_and_genebody_beds(
    bed_nonoverlap,
    tsr_final,
    gb_final,
    config$tsr_parameters$upstream_extension,
    config$tsr_parameters$downstream_extension,
    gb_extension
  )

  # Step 8: Compute signals over TSR and gene body
  log_info("\n### Step 8: Computing signals over TSR and gene body ###")
  tsr_signal <- file.path(signal_dir, paste0(prefix, "_tsr_signal.txt"))
  gb_signal <- file.path(signal_dir, paste0(prefix, "_genebody_signal.txt"))

  compute_signal_over_bed(
    config$input$sample_bigwig,
    tsr_final,
    tsr_signal,
    config$tools$bigwig_average_over_bed
  )

  compute_signal_over_bed(
    config$input$sample_bigwig,
    gb_final,
    gb_signal,
    config$tools$bigwig_average_over_bed
  )

  # Step 9: Compute Travel Ratio
  log_info("\n### Step 9: Computing Travel Ratio ###")
  tr_output <- file.path(output_dir, paste0(prefix, "_travel_ratio.txt"))
  tr_data <- compute_travel_ratio(tsr_signal, gb_signal, tr_output)

  # Step 10: Plot Travel Ratio distribution
  if (config$plotting$generate_plots) {
    log_info("\n### Step 10: Generating Travel Ratio plots ###")
    tr_plot <- file.path(plot_dir, paste0(prefix, "_travel_ratio_cdf.",
                                         config$plotting$plot_format))
    plot_travel_ratio_cdf(tr_output, tr_plot, config)
  }

  # Cleanup intermediate files if requested
  if (!config$output$keep_intermediate) {
    log_info("\n### Cleaning up intermediate files ###")
    intermediate_files <- c(tsr_for_selection, selection_signal,
                           unique_transcripts, bed_unique)
    for (file in intermediate_files) {
      if (file.exists(file)) {
        file.remove(file)
        log_debug("Removed: {file}")
      }
    }
    log_success("Intermediate files removed")
  }

  # Final summary
  end_time <- Sys.time()
  elapsed_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  log_info("\n" %R% strrep("=", 72))
  log_success("Pipeline completed successfully!")
  log_info("Total runtime: {round(elapsed_time, 1)} seconds")
  log_info("Output directory: {output_dir}")
  log_info("Travel Ratio file: {tr_output}")
  log_info(strrep("=", 72))

  # Return results
  return(list(
    travel_ratio_file = tr_output,
    output_directory = output_dir,
    final_bed = bed_nonoverlap,
    n_genes = nrow(tr_data),
    runtime_seconds = elapsed_time
  ))
}

#===============================================================================
# EXECUTION
#===============================================================================

# Run pipeline if script is executed directly
if (!interactive()) {

  # Handle config file override
  config_file <- opt$config

  if (!is.null(opt$output)) {
    # Load config and override output directory
    config <- yaml::read_yaml(config_file)
    config$output$output_dir <- opt$output

    # Write temporary config
    temp_config <- tempfile(fileext = ".yaml")
    yaml::write_yaml(config, temp_config)
    config_file <- temp_config
  }

  # Set verbosity
  if (opt$verbose) {
    log_threshold(DEBUG)
  }

  # Run pipeline
  tryCatch({
    results <- run_travel_ratio_pipeline(config_file)

    cat("\n")
    cat("Travel Ratio computation completed successfully!\n")
    cat("Results saved to:", results$travel_ratio_file, "\n")
    cat("Number of genes analyzed:", results$n_genes, "\n")

  }, error = function(e) {
    log_error("Pipeline failed with error: {e$message}")
    cat("\nERROR: Pipeline failed. Check log file for details.\n", file = stderr())
    quit(status = 1)
  })

  # Clean up temp config if created
  if (exists("temp_config") && file.exists(temp_config)) {
    file.remove(temp_config)
  }
}
