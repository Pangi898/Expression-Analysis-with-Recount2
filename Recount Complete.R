## Load libraries
{
library(data.table)
library(recount)
library(dplyr)
library(tibble) 
library(writexl)
library(ggplot2)
library(rstatix)
library(gtools)
library(purrr)
library(rlang)
}
  
## Read ensembl coordinates file
coordinates <- fread("D:/University/Bachelor/Διπλωματική Εργασία/Coordinates/Gene Coordinates.txt")
#view (coordinates)
## Rename columns and chromosomes in the coordinates data frame
coordinates_renamed <- coordinates %>%
  rename(
    ENSG_ID = `Gene stable ID version`,
    chrom = `Chromosome/scaffold name`,
    geneStart = `Gene start (bp)`,
    geneEnd = `Gene end (bp)`,
    gene_name = `Gene name`
  ) %>%
  mutate(chrom = case_when(
    chrom == "1" ~ "01",
    chrom == "2" ~ "02",
    chrom == "3" ~ "03",
    chrom == "4" ~ "04",
    chrom == "5" ~ "05",
    chrom == "6" ~ "06",
    chrom == "7" ~ "07",
    chrom == "8" ~ "08",
    chrom == "9" ~ "09",
    TRUE ~ chrom
  ))
{
#view(coordinates_renamed)
##import chromosome lengths 
chrom_lengths<- fread("D:/University/Bachelor/Διπλωματική Εργασία/Coordinates/hg38 - chrom lengths.genome")
#view(Chrom_lengths)
chrom_lengths_renamed <- chrom_lengths %>%
  rename(
    chrom = V1,
    length = V2
  ) %>%
  mutate(chrom = case_when(
    chrom == "chr1" ~ "01",
    chrom == "chr2" ~ "02",
    chrom == "chr3" ~ "03",
    chrom == "chr4" ~ "04",
    chrom == "chr5" ~ "05",
    chrom == "chr6" ~ "06",
    chrom == "chr7" ~ "07",
    chrom == "chr8" ~ "08",
    chrom == "chr9" ~ "09",
    chrom == "chr10" ~ "10",
    chrom == "chr11" ~ "11",
    chrom == "chr12" ~ "12",
    chrom == "chr13" ~ "13",
    chrom == "chr14" ~ "14",
    chrom == "chr15" ~ "15",
    chrom == "chr16" ~ "16",
    chrom == "chr17" ~ "17",
    chrom == "chr18" ~ "18",
    chrom == "chr19" ~ "19",
    chrom == "chr20" ~ "20",
    chrom == "chr21" ~ "21",
    chrom == "chr22" ~ "22",
    chrom == "chrX" ~ "X",
    chrom == "chrY" ~ "Y",
    TRUE ~ chrom
  ))
#view(chrom_lengths_renamed)
}
## Find the single - cell project IDs by searching abstracts of studies
# abstract_search("single cell")

## Categorize the studies
human_stem_cells_studies <- c("SRP011546", "SRP033135", "SRP034712", "SRP040145", "SRP050992", "SRP053004",
                              "SRP056957", "SRP056989", "SRP058773", "SRP058977", "SRP059379", "SRP063754",
                              "SRP066834", "SRP066994")

human_cancer_cells_studies <- c("DRP001358", "DRP002586", "ERP000959", "ERP006662", "SRP013299", "SRP014842",
                                "SRP018838", "SRP022259", "SRP022260", "SRP030617", "SRP040288", "SRP042161",
                                "SRP048222", "SRP057508", "SRP058046", "SRP060416", "SRP063840")

human_differentiated_cells_studies <- c("SRP048971", "SRP055153", "SRP057196")

## Create empty lists to store the expression counts data frames for each category
stem_cell_expression_data_list <- list()
cancer_cell_expression_data_list <- list()
differentiated_cell_expression_data_list <- list()

##Function to process studies and return transcript-aggrevated expression data 
process_studies <- function(study_vector) {
  expression_data_list <- list()
  
  for (study_id in study_vector) {
    
    ## Define the file path for the current study dynamically
    file_path <- file.path(study_id, "rse_gene.Rdata")
    
    ## Download the data if it is not there
    if (!file.exists(file_path)) {
      download_study(study_id, type = "rse-gene")
    }
    
    ## Check if the file exists
    if (file.exists(file_path)) {
      
      ## Load the data
      load(file_path, verbose = TRUE)
      
      ## Check if the rse_gene object exists
      if (exists("rse_gene")) {
        
        ## Scale the counts to RPM using the scale_counts function from recount
        rse_gene_scaled <- scale_counts(rse_gene)
        
        ## Extract the scaled expression matrix
        expression_matrix_scaled <- as.data.frame(assay(rse_gene_scaled))
        
        ## Convert rownames (ESNG IDs) to a column named 'ENSG_ID'
        expression_matrix_scaled <- expression_matrix_scaled %>%
          rownames_to_column(var = "ENSG_ID")
        
        ## Store the resulting expression matrix in the list with study_id as the key
        expression_data_list[[study_id]] <- expression_matrix_scaled
        
        ## Print a message indicating successful processing of the current study
        message(paste("Processed study:", study_id))
        
      } else {
        warning(paste("rse_gene object not found for", study_id))
      }
      
    } else {
      warning(paste("Data for", study_id, "could not be found."))
    }
  }
  
  return(expression_data_list)
}

## Process each vector separately
stem_cell_expression_data_list <- process_studies(human_stem_cells_studies)
cancer_cell_expression_data_list <- process_studies(human_cancer_cells_studies)
differentiated_cell_expression_data_list <- process_studies(human_differentiated_cells_studies)

## Checkpoint 1: 
#view(stem_cell_expression_data_list$SRP066994)
#view(cancer_cell_expression_data_list$DRP001358)
#view(differentiated_cell_expression_data_list$SRP048971)

##Function to filter genes with 0 expression in at least 10 samples
filter_genes <- function(expression_data) {
  expression_data_filtered <- expression_data %>%
    ## Check how many zeros each gene has
    rowwise() %>%
    ## Count zeros per gene across all samples
    mutate(zero_count = sum(across(-ENSG_ID, ~ . == 0))) %>%
    ## Filter genes with at least 10 zeros
    filter(zero_count < 10) %>%
    ## Drop the helper zero_count column
    select(-zero_count)
  
  return(expression_data_filtered)
}

## Apply the filtering function to every data frame in each list
stem_cell_expression_data_list <- lapply(stem_cell_expression_data_list, filter_genes)
cancer_cell_expression_data_list <- lapply(cancer_cell_expression_data_list, filter_genes)
differentiated_cell_expression_data_list <- lapply(differentiated_cell_expression_data_list, filter_genes)

## Checkpoint 2:
#view(stem_cell_expression_data_list$SRP066994)

## Initialize empty lists to store the average expression for each study in each category
average_stem_cell_expression_list <- list()
average_cancer_cell_expression_list <- list()
average_differentiated_cell_expression_list <- list()

## Function to calculate average expression for each gene in each study
calculate_average_expression <- function(expression_data_list) {
  average_expression_list <- list()
  
  ## Loop through each study in expression_data_list
  for (study_id in names(expression_data_list)) {
    
    ## Extract the study data
    study_data <- expression_data_list[[study_id]]
    
    ## Calculate the average expression for each gene in the current study
    average_expression_df <- study_data %>%
      rowwise() %>%  ## Apply function to each row (gene)
      summarise(
        ENSG_ID = ENSG_ID,  ## Keep the ENSG_ID
        ## Calculate the mean excluding 0 values, convert to percentage
        average_expression_counts = mean(c_across(-ENSG_ID)[c_across(-ENSG_ID) > 0], na.rm = TRUE) * 100
      ) %>%
      ## Filter out genes where average expression is less than 1%
      filter(average_expression_counts >= 1)
    
    ## Store the result in the list with study ID as the key
    average_expression_list[[study_id]] <- average_expression_df
  }
  
  return(average_expression_list)
}

## Apply the function to each of the three expression lists
average_stem_cell_expression_list <- calculate_average_expression(stem_cell_expression_data_list)
average_cancer_cell_expression_list <- calculate_average_expression(cancer_cell_expression_data_list)
average_differentiated_cell_expression_list <- calculate_average_expression(differentiated_cell_expression_data_list)

## Checkpoint 3:
#view(average_stem_cell_expression_list$SRP066994)

## Studies with the most genes from each category: SRP066994, SRP048222, SRP048971

## Function to merge average expression data for a single study with coordinates
merge_single_with_coordinates <- function(study_id, average_expression_list, coordinates_df) {
  if (!study_id %in% names(average_expression_list)) {
    stop(paste("Study ID", study_id, "not found in the expression list."))
  }
  
  ## Merge the expression data with coordinates
  merged_df <- average_expression_list[[study_id]] %>%
    inner_join(coordinates_df, by = c("ENSG_ID" = "ENSG_ID")) %>%
    arrange(chrom, geneStart)  ## Reorder by chromosome and gene start
  
  return(merged_df)
}

## Call the merge_single_with_coordinates function for SRP066994, SRP048222, SRP048971
SRP066994_merged_df <- merge_single_with_coordinates("SRP066994", average_stem_cell_expression_list, coordinates_renamed)
SRP048222_merged_df <- merge_single_with_coordinates("SRP048222", average_cancer_cell_expression_list, coordinates_renamed)
SRP048971_merged_df <- merge_single_with_coordinates("SRP048971", average_differentiated_cell_expression_list, coordinates_renamed)
## Checkpoint 4:
#view(SRP066994_merged_df)
#view(SRP048222_merged_df)
#view(SRP048971_merged_df)
## Function to divide the avg score with the gene length (geneEnd - geneStart) for a single DataFrame
divide_expression_with_distance_function <- function(data) {
  ## Ensure required columns exist
  required_columns <- c("average_expression_counts", "geneStart", "geneEnd")
  
  if (!all(required_columns %in% colnames(data))) {
    stop("One or more required columns are missing from the data frame.")
  }
  
  ## Perform the calculation and update the column
  modified_data <- data %>%
    mutate(average_expression_counts = average_expression_counts / (geneEnd - geneStart))
  
  ## Return the modified data frame
  return(modified_data)
}

## Call the divide_expression_with_distance_function function for SRP066994, SRP048222, SRP048971
SRP066994_length_standarized_df <- divide_expression_with_distance_function(SRP066994_merged_df)
SRP048222_length_standarized_df <- divide_expression_with_distance_function(SRP048222_merged_df)
SRP048971_length_standarized_df <- divide_expression_with_distance_function(SRP048971_merged_df)
## Checkpoint 5:
#view(SRP066994_length_standarized_df)
#view(SRP048222_length_standarized_df)
#view(SRP048971_length_standarized_df)
# Function to split a dataframe by the 'chrom' column into a list of data frames
split_by_chrom <- function(df) {
  # Check if the 'chrom' column exists in the dataframe
  if (!"chrom" %in% colnames(df)) {
    stop("The dataframe must contain a 'chrom' column.")
  }
  
  # Split the dataframe based on 'chrom' column into a list of data frames
  chromosome_split_list <- df %>%
    group_by(chrom) %>%
    group_split()
  
  # Assign names to each data frame in the list based on the 'chrom' values
  names(chromosome_split_list) <- unique(df$chrom)
  
  # Return the list of data frames
  return(chromosome_split_list)
}

## Call the split_by_chrom function for SRP066994, SRP048222, SRP048971
SRP066994_chromosome_split_list <- split_by_chrom(SRP066994_length_standarized_df)
SRP048222_chromosome_split_list <- split_by_chrom(SRP048222_length_standarized_df)
SRP048971_chromosome_split_list <- split_by_chrom(SRP048971_length_standarized_df)

## Checkpoint 6:
#view(SRP066994_chromosome_split_list$`15`)
#view(SRP048222_chromosome_split_list$'15')
#view(SRP048971_chromosome_split_list$'15')

## Function to calculate expression differences with gene midpoint
calculate_expression_differences <- function(
    chromosome_split_list, 
    distance_threshold = NULL, 
    num_genes = NULL,
    seed = NULL
) {
  if (!is.null(seed)) set.seed(seed)
  
  expression_differences_list <- list()
  
  for (chrom in names(chromosome_split_list)) {
    chrom_df <- chromosome_split_list[[chrom]]
    
    midpoints <- c()
    neighborhood_differences <- c()
    control_differences <- c()
    
    if (nrow(chrom_df) > 1) {
      expx <- chrom_df$average_expression_counts
      gene_start <- chrom_df$geneStart
      gene_end <- chrom_df$geneEnd
      gene_positions <- (gene_start + gene_end) / 2
      
      for (i in seq_along(expx)) {
        if (!is.null(distance_threshold)) {
          distances <- abs(gene_positions - gene_positions[i])
          neighbors <- which(distances > 0 & distances <= distance_threshold)
        } else if (!is.null(num_genes)) {
          half_window <- floor((num_genes - 1) / 2)
          neighbors <- seq(max(1, i - half_window), min(length(expx), i + half_window))
          neighbors <- setdiff(neighbors, i)
          if (length(neighbors) < (num_genes - 1)) next
        } else {
          stop("Either distance_threshold or num_genes must be specified.")
        }
        
        if (length(neighbors) == 0) next
        
        surrounding_mean <- mean(expx[neighbors], na.rm = TRUE)
        exp_diff <- abs(expx[i] - surrounding_mean)
        
        available_indices <- setdiff(seq_along(expx), c(i, neighbors))
        if (!is.null(distance_threshold)) {
          control_neighbors <- which(abs(gene_positions - gene_positions[i]) > distance_threshold)
        } else if (!is.null(num_genes) && length(available_indices) >= length(neighbors)) {
          control_neighbors <- sample(available_indices, length(neighbors))
        } else {
          control_neighbors <- NULL
        }
        
        if (is.null(control_neighbors) || length(control_neighbors) == 0) next
        
        control_mean <- mean(expx[control_neighbors], na.rm = TRUE)
        control_diff <- abs(expx[i] - control_mean)
        
        if (!is.na(exp_diff) && !is.na(control_diff)) {
          midpoints <- c(midpoints, gene_positions[i])
          neighborhood_differences <- c(neighborhood_differences, exp_diff)
          control_differences <- c(control_differences, control_diff)
        }
      }
    }
    
    expression_differences_list[[chrom]] <- data.frame(
      midpoint = midpoints,
      neighborhood_differences = neighborhood_differences,
      control_differences = control_differences
    )
  }
  
  return(expression_differences_list)
}

##Apply the calculate_expression_differences for individual distances or gene_numbers
#SRP066994_1000000_expr_diffs_list <- calculate_expression_differences(SRP066994_chromosome_split_list, distance_threshold = 1000000, seed = 123)
#SRP066994_3_expr_diffs_list <- calculate_expression_differences(SRP066994_chromosome_split_list, distance_threshold = NULL, num_genes = 3, seed = 123)
SRP048971_1000000_expr_diffs_list <- calculate_expression_differences(SRP048971_chromosome_split_list, distance_threshold = 1000000, seed = 123)
SRP048971_3_expr_diffs_list <- calculate_expression_differences(SRP048971_chromosome_split_list, distance_threshold = NULL, num_genes = 3, seed = 123)

##Checkpoint 7
#view(SRP066994_1000000_expr_diffs_list$`01`)
#view(SRP066994_3_expr_diffs_list$`01`)
view(SRP048971_1000000_expr_diffs_list$`01`)
view(SRP048971_3_expr_diffs_list$`01`)

## Function to calculate ratios per midpoint
calculate_ratios_per_midpoint <- function(expression_differences_list) {
  ## Initialize result list
  ratios_list <- list()
  
  ## Loop through each chromosome
  for (chrom in names(expression_differences_list)) {
    df <- expression_differences_list[[chrom]]
    
    ## Only proceed if both difference columns exist and have rows
    if (
      "neighborhood_differences" %in% colnames(df) &&
      "control_differences" %in% colnames(df) &&
      nrow(df) > 0
    ) {
      ## Compute the ratio; use NA for division by zero
      df$ratios <- with(df, ifelse(control_differences != 0,
                                   neighborhood_differences / control_differences,
                                   NA))
    } else {
      df$ratios <- numeric(0)  # empty column if the data frame is empty
    }
    
    ## Store in result list
    ratios_list[[chrom]] <- df
  }
  
  return(ratios_list)
}

##Apply the calculate_ratios_per_midpoint for individual distances or gene_numbers
#SRP066994_1000000_midpoint_ratios <- calculate_ratios_per_midpoint(SRP066994_1000000_expr_diffs_list)
#SRP066994_3_midpoint_ratios <- calculate_ratios_per_midpoint(SRP066994_3_expr_diffs_list)
SRP048971_1000000_midpoint_ratios <- calculate_ratios_per_midpoint(SRP048971_1000000_expr_diffs_list)
SRP048971_3_midpoint_ratios <- calculate_ratios_per_midpoint(SRP048971_3_expr_diffs_list)

##Chechpoint 8
#view(SRP066994_1000000_midpoint_ratios$`01`)
#view(SRP066994_3_midpoint_ratios$`01`)
view(SRP048971_1000000_midpoint_ratios$`01`)
view(SRP048971_3_midpoint_ratios$`01`)

## function to plot ratios across midpoints, faceted by chromosome (change title depending on the parameter)
plot_ratios_across_midpoints_facet <- function(ratios_per_midpoint_data) {
  combined_df <- do.call(rbind, lapply(names(ratios_per_midpoint_data), function(chrom) {
    df <- ratios_per_midpoint_data[[chrom]]
    if (!is.null(df$ratios) && length(df$ratios) > 0) {
      df$midpoint_index <- seq_len(nrow(df))
      df$chromosome <- chrom  # Add chromosome as a new column
      return(df)
    } else {
      return(NULL)
    }
  }))
  
  if (is.null(combined_df)) return(NULL)  # Return early if no valid data
  
  p <- ggplot(combined_df, aes(x = midpoint_index, y = log10(ratios))) +
    geom_col(fill = "steelblue") +
    facet_wrap(~ chromosome, scales = "free_x") +
    labs(
      title = "log10(Ratio) Across Midpoints by Chromosome - 3",
      x = "Midpoint Index",
      y = "log10(Ratio)"
    ) +
    theme_minimal()
  
  return(p)
}

##Apply the plot_ratios_across_midpoints function 
#SRP066994_1000000_midpoint_plots <- plot_ratios_across_midpoints_facet(SRP066994_1000000_midpoint_ratios)
#SRP066994_3_midpoint_plots <- plot_ratios_across_midpoints_facet(SRP066994_3_midpoint_ratios)
#SRP048971_1000000_midpoint_plots <- plot_ratios_across_midpoints_facet(SRP048971_1000000_midpoint_ratios)
SRP048971_3_midpoint_plots <- plot_ratios_across_midpoints_facet(SRP048971_3_midpoint_ratios)

##Checkpoint 9
#print(SRP066994_1000000_midpoint_plots)
#print(SRP066994_3_midpoint_plots)
#print(SRP048971_1000000_midpoint_plots)
print(SRP048971_3_midpoint_plots)

##function to calculate means ratios and p-values
calculate_means_ratios_and_p_values <- function(expression_differences_list) {
  results <- data.frame(
    chromosome = character(),
    mean_ratio = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (chrom in names(expression_differences_list)) {
    df <- expression_differences_list[[chrom]]
    
    if (!is.null(df) &&
        nrow(df) > 0 &&
        all(c("neighborhood_differences", "control_differences") %in% names(df))) {
      
      neighborhood <- df$neighborhood_differences
      control <- df$control_differences
      
      # Ensure both columns have values
      if (length(neighborhood) > 0 && length(control) > 0) {
        neighborhood_mean <- mean(neighborhood, na.rm = TRUE)
        control_mean <- mean(control, na.rm = TRUE)
        
        # Skip if control mean is 0 to avoid division by zero
        if (control_mean == 0 || is.na(control_mean)) next
        
        ratio <- neighborhood_mean / control_mean
        
        # Perform Wilcoxon rank-sum test
        p_val <- tryCatch({
          wilcox.test(neighborhood, control)$p.value
        }, error = function(e) {
          NA
        })
        
        results <- rbind(
          results,
          data.frame(
            chromosome = chrom,
            mean_ratio = ratio,
            p_value = p_val,
            stringsAsFactors = FALSE
          )
        )
      }
    }
  }
  
  return(results)
}

##Apply the calculate_means_ratios_and_p_values for individual distances
SRP66994_10000000_means_ratios_and_p_values <- calculate_means_ratios_and_p_values(SRP066994_10000000_expr_diffs_list)
SRP66994_1000000_means_ratios_and_p_values <- calculate_means_ratios_and_p_values(SRP066994_1000000_expr_diffs_list)
SRP66994_100000_means_ratios_and_p_values <- calculate_means_ratios_and_p_values(SRP066994_100000_expr_diffs_list)

##Checkpoint 10
#view(SRP66994_1000000_means_ratios_and_p_values)
#view(SRP66994_100000_means_ratios_and_p_values_arranged)

##plot mean ratios for individual distances (change title based on distance)
plot_mean_ratios <- function(df) {
  p1 <- ggplot(df, aes(x = chromosome, y = mean_ratio)) +
    geom_point(color = "blue") +
    geom_hline(yintercept = 1, color = "red", linetype = "dotted", linewidth = 1) +
    labs(title = "Mean Ratio vs Chromosome, ...bp",
         x = "Chromosome",
         y = "Mean Ratio") +
    theme_minimal()
  
  return(p1)
}

##Apply the plot_mean_ratios_and_p_values for individual distances for all chromosomes
#SRP66994_10000000_mean_ratio_plot <- plot_mean_ratios(SRP66994_10000000_means_ratios_and_p_values)
#SRP66994_1000000_mean_ratio_plot <- plot_mean_ratios(SRP66994_1000000_means_ratios_and_p_values)
#SRP66994_100000_mean_ratio_plot <- plot_mean_ratios(SRP66994_100000_means_ratios_and_p_values)

##Checkpoint 11
#print(SRP66994_10000000_mean_ratio_plot)
#print(SRP66994_1000000_mean_ratio_plot)
#print(SRP66994_100000_mean_ratio_plot)

##function to merge the means ratios and p-values dfs with chrom_lengths_renamed
merge_with_chrom_lengths <- function(df) {
  # Ensure Chrom_lengths_renamed exists in the environment
  if (!exists("Chrom_lengths_renamed")) {
    stop("Chrom_lengths_renamed does not exist in the environment.")
  }
  
  # Perform the merge using left_join
  merged_df <- df %>%
    left_join(Chrom_lengths_renamed, by = c("chromosome" = "chrom"))
  
  return(merged_df)
}

#apply the merge_with_chrom_lengths function for individual distances or gene numbers
SRP66994_1000000_means_ratios_and_p_values_merged <- merge_with_chrom_lengths(SRP66994_1000000_means_ratios_and_p_values)

##Checkpoint 12
#view(SRP66994_1000000_means_ratios_and_p_values_merged)

#function to plot mean ratios across lengths
plot_means_ratios_across_lengths <- function(df) {
  # Basic checks
  required_cols <- c("length", "mean_ratio", "chromosome")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }
  
  # Create the plot
  p <- ggplot(df, aes(x = length, y = mean_ratio, color = chromosome)) +
    geom_point(size = 2, alpha = 0.7) +
    labs(
      title = "Mean Ratios Across Chromosome Lengths",
      x = "Chromosome Length",
      y = "Mean Ratio"
    ) +
    theme_minimal() +
    theme(
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    )
  
  return(p)
}

##apply the plot_means_ratios_across_lengths
SRP66994_1000000_means_ratios_across_length_plot <- plot_means_ratios_across_lengths(SRP66994_1000000_means_ratios_and_p_values_merged)
##Checkpoint 13
print(SRP66994_1000000_means_ratios_across_length_plot)

##DIFFERING DISTANCES & GENE NUMBERS
## Apply the calculate_expression_differences function for SRP066994, SRP048222, SRP048971 for differing distances
## Create a vector with distances for iteration
distances_vector <- c(100000, 200000, 300000, 400000, 500000, 600000, 700000, 800000, 900000, 1000000, 2000000, 3000000, 4000000, 5000000, 6000000, 7000000, 8000000, 9000000, 10000000)
##SRP066994
# Initialize an empty list to store results
SRP66994_differing_distances_expr_diffs_nested_list <- list()
# Loop through each distance and apply the function
for (dist in distances_vector) {
  message("Processing distance: ", dist)
  
  result <- calculate_expression_differences(
    chromosome_split_list = SRP066994_chromosome_split_list,
    distance_threshold = dist, 
    seed = 123
  )
  
  SRP66994_differing_distances_expr_diffs_nested_list[[as.character(dist)]] <- result
}

##SRP048222
# Initialize an empty list to store results
SRP048222_differing_distances_expr_diffs_nested_list <- list()
# Loop through each distance and apply the function
for (dist in distances_vector) {
  message("Processing distance: ", dist)
  
  result <- calculate_expression_differences(
    chromosome_split_list = SRP048222_chromosome_split_list,
    distance_threshold = dist, 
    seed = 123
  )
  
  SRP048222_differing_distances_expr_diffs_nested_list[[as.character(dist)]] <- result
}
##SRP048971
# Initialize an empty list to store results
SRP048971_differing_distances_expr_diffs_nested_list <- list()
# Loop through each distance and apply the function
for (dist in distances_vector) {
  message("Processing distance: ", dist)
  
  result <- calculate_expression_differences(
    chromosome_split_list = SRP048971_chromosome_split_list,
    distance_threshold = dist, 
    seed = 123
  )
  
  SRP048971_differing_distances_expr_diffs_nested_list[[as.character(dist)]] <- result
}

## Checkpoint 14: 
#view(SRP66994_differing_distances_expr_diffs_nested_list$`1e+06`$`01`)
#view(SRP048222_differing_distances_expr_diffs_nested_list$`1e+06`$`01`)
#view(SRP048971_differing_distances_expr_diffs_nested_list$`1e+06`$`01`)

## Apply the calculate_means_ratios_and_p_values function for SRP066994, SRP048222, SRP048971 for differing distances
##SRP066994
# Initialize an empty list to store results
SRP066994_differing_distances_means_ratios_and_pvalues_list <- list()

# Loop through each sublist (each distance)
for (dist in names(SRP66994_differing_distances_expr_diffs_nested_list)) {
  message("Calculating statistics for distance: ", dist)
  
  result_df <- calculate_means_ratios_and_p_values(SRP66994_differing_distances_expr_diffs_nested_list[[dist]])
  
  SRP066994_differing_distances_means_ratios_and_pvalues_list[[dist]] <- result_df
}

##SRP048222
# Initialize an empty list to store results
SRP048222_differing_distances_means_ratios_and_pvalues_list <- list()

# Loop through each sublist (each distance)
for (dist in names(SRP048222_differing_distances_expr_diffs_nested_list)) {
  message("Calculating statistics for distance: ", dist)
  
  result_df <- calculate_means_ratios_and_p_values(SRP048222_differing_distances_expr_diffs_nested_list[[dist]])
  
  SRP048222_differing_distances_means_ratios_and_pvalues_list[[dist]] <- result_df
}
##SRP048971
# Initialize an empty list to store results
SRP048971_differing_distances_means_ratios_and_pvalues_list <- list()

# Loop through each sublist (each distance)
for (dist in names(SRP048971_differing_distances_expr_diffs_nested_list)) {
  message("Calculating statistics for distance: ", dist)
  
  result_df <- calculate_means_ratios_and_p_values(SRP048971_differing_distances_expr_diffs_nested_list[[dist]])
  
  SRP048971_differing_distances_means_ratios_and_pvalues_list[[dist]] <- result_df
}

## Checkpoint 15:
#view(SRP066994_differing_distances_means_ratios_and_pvalues_list$`1e+06`)
#view(SRP048222_differing_distances_means_ratios_and_pvalues_list$`1e+06`)
#view(SRP048971_differing_distances_means_ratios_and_pvalues_list$`1e+06`)

## Function to calculate average statistics per distance with Standard Error
calculate_average_statistics_per_distance <- function(data_list) {
  # Loop through each distance (name of the list element)
  result_list <- lapply(names(data_list), function(dist_name) {
    df <- data_list[[dist_name]]
    
    # Calculate sample sizes (excluding NAs)
    n_mean_ratio <- sum(!is.na(df$mean_ratio))
    n_p_value    <- sum(!is.na(df$p_value))
    
    # Calculate average and standard error for mean_ratio and p_value
    avg_mean_ratio <- mean(df$mean_ratio, na.rm = TRUE)
    se_mean_ratio  <- sd(df$mean_ratio, na.rm = TRUE) / sqrt(n_mean_ratio)
    
    avg_p_value <- mean(df$p_value, na.rm = TRUE)
    se_p_value  <- sd(df$p_value, na.rm = TRUE) / sqrt(n_p_value)
    
    # Create a row with all statistics
    data.frame(
      distance = as.numeric(dist_name),
      average_mean_ratio = avg_mean_ratio,
      SE_mean_ratio = se_mean_ratio,
      average_p_value = avg_p_value,
      SE_p_value = se_p_value
    )
  })
  
  # Combine all rows into a single dataframe
  result_df <- do.call(rbind, result_list)
  rownames(result_df) <- NULL
  return(result_df)
}

##Apply the calculate_average_statistics_per_distance function for SRP066994, SRP048222, SRP048971
SRP066994_differing_distances_statistics_averages <- calculate_average_statistics_per_distance(SRP066994_differing_distances_means_ratios_and_pvalues_list)
SRP048222_differing_distances_statistics_averages <- calculate_average_statistics_per_distance(SRP048222_differing_distances_means_ratios_and_pvalues_list)
SRP048971_differing_distances_statistics_averages <- calculate_average_statistics_per_distance(SRP048971_differing_distances_means_ratios_and_pvalues_list)

##Checkpoint 16
#view(SRP066994_differing_distances_statistics_averages)
#view(SRP048222_differing_distances_statistics_averages)
#view(SRP048971_differing_distances_statistics_averages)

# Function to plot average mean ratios and p-values across distances with SE bands
plot_average_mean_ratios_and_p_values_distances <- function(df) {
  # Plot 1: average_mean_ratio vs distance with SE bands
  p1 <- ggplot(df, aes(x = distance, y = average_mean_ratio)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    geom_ribbon(aes(ymin = average_mean_ratio - SE_mean_ratio,
                    ymax = average_mean_ratio + SE_mean_ratio),
                alpha = 0.2, fill = "blue") +
    labs(title = "Average Mean Ratio vs Distance",
         x = "Distance",
         y = "Average Mean Ratio") +
    theme_minimal()
  
  # Plot 2: average_p_value vs distance with SE bands
  p2 <- ggplot(df, aes(x = distance, y = average_p_value)) +
    geom_line(color = "red") +
    geom_point(color = "red") +
    geom_ribbon(aes(ymin = average_p_value - SE_p_value,
                    ymax = average_p_value + SE_p_value),
                alpha = 0.2, fill = "red") +
    labs(title = "Average p-Value vs Distance",
         x = "Distance",
         y = "Average p-Value") +
    theme_minimal()
  
  return(list(mean_ratio_plot = p1, p_value_plot = p2))
}

##Apply the plot_average_mean_ratios_and_p_values for SRP066994, SRP048222, SRP048971
SRP066994_distances_plots_list <- plot_average_mean_ratios_and_p_values_distances(SRP066994_differing_distances_statistics_averages)
SRP048222_distances_plots_list <- plot_average_mean_ratios_and_p_values_distances(SRP048222_differing_distances_statistics_averages)
SRP048971_distances_plots_list <- plot_average_mean_ratios_and_p_values_distances(SRP048971_differing_distances_statistics_averages)

##Checkpoint 17
#print(SRP066994_distances_plots_list$mean_ratio_plot)
#print(SRP066994_distances_plots_list$p_value_plot)
#print(SRP048222_distances_plots_list$mean_ratio_plot)
#print(SRP048222_distances_plots_list$p_value_plot)
#print(SRP048971_distances_plots_list$mean_ratio_plot)
#print(SRP048971_distances_plots_list$p_value_plot)

## Apply the calculate_expression_differences function for SRP066994, SRP048222, SRP048971 for differing gene numbers
## Create vectors with gene numbers for iteration
gene_numbers_vector <- c(3, 13, 23, 33, 43, 53, 63, 73, 83, 93, 103)

##SRP066994
# Initialize an empty list to store results
SRP66994_differing_num_genes_expr_diffs_nested_list <- list()
# Loop through each distance and apply the function
for (num in gene_numbers_vector) {
  message("Processing gene number: ", num)
  
  result <- calculate_expression_differences(
    chromosome_split_list = SRP066994_chromosome_split_list,
    distance_threshold = NULL,
    num_genes = num, 
    seed = 123
  )
  
  SRP66994_differing_num_genes_expr_diffs_nested_list[[as.character(num)]] <- result
}

##SRP048222
# Initialize an empty list to store results
SRP048222_differing_num_genes_expr_diffs_nested_list <- list()
# Loop through each distance and apply the function
for (num in gene_numbers_vector) {
  message("Processing gene number: ", num)
  
  result <- calculate_expression_differences(
    chromosome_split_list = SRP048222_chromosome_split_list,
    distance_threshold = NULL,
    num_genes = num, 
    seed = 123
  )
  
  SRP048222_differing_num_genes_expr_diffs_nested_list[[as.character(num)]] <- result
}

##SRP048971
# Initialize an empty list to store results
SRP048971_differing_num_genes_expr_diffs_nested_list <- list()
# Loop through each distance and apply the function
for (num in gene_numbers_vector) {
  message("Processing gene number: ", num)
  
  result <- calculate_expression_differences(
    chromosome_split_list = SRP048971_chromosome_split_list,
    distance_threshold = NULL,
    num_genes = num, 
    seed = 123
  )
  
  SRP048971_differing_num_genes_expr_diffs_nested_list[[as.character(num)]] <- result
}

## Checkpoint 18: 
#view(SRP66994_differing_num_genes_expr_diffs_nested_list$`23`$`02`)
#view(SRP048222_differing_num_genes_expr_diffs_nested_list$`23`$`02`)
#view(SRP048971_differing_num_genes_expr_diffs_nested_list$`23`$`02`)

## Apply the calculate_means_ratios_and_p_values function for SRP066994, SRP048222, SRP048971 for gene numbers
##SRP066994 
# Initialize an empty list to store results
SRP066994_differing_num_genes_means_ratios_and_pvalues_list<- list()
# Loop through each sublist 
for (num in names(SRP66994_differing_num_genes_expr_diffs_nested_list)) {
  message("Calculating statistics for gene number: ", num)
  
  result_df <- calculate_means_ratios_and_p_values(SRP66994_differing_num_genes_expr_diffs_nested_list[[num]])
  
  SRP066994_differing_num_genes_means_ratios_and_pvalues_list[[num]] <- result_df
}

##SRP048222
# Initialize an empty list to store results
SRP048222_differing_num_genes_means_ratios_and_pvalues_list<- list()
# Loop through each sublist 
for (num in names(SRP048222_differing_num_genes_expr_diffs_nested_list)) {
  message("Calculating statistics for gene number: ", num)
  
  result_df <- calculate_means_ratios_and_p_values(SRP048222_differing_num_genes_expr_diffs_nested_list[[num]])
  
  SRP048222_differing_num_genes_means_ratios_and_pvalues_list[[num]] <- result_df
}
##SRP048971
# Initialize an empty list to store results
SRP048971_differing_num_genes_means_ratios_and_pvalues_list<- list()
# Loop through each sublist 
for (num in names(SRP048971_differing_num_genes_expr_diffs_nested_list)) {
  message("Calculating statistics for gene number: ", num)
  
  result_df <- calculate_means_ratios_and_p_values(SRP048971_differing_num_genes_expr_diffs_nested_list[[num]])
  
  SRP048971_differing_num_genes_means_ratios_and_pvalues_list[[num]] <- result_df
}
## Checkpoint 19:
#view(SRP066994_differing_num_genes_means_ratios_and_pvalues_list$`23`)
#view(SRP048222_differing_num_genes_means_ratios_and_pvalues_list$`23`)
#view(SRP048971_differing_num_genes_means_ratios_and_pvalues_list$`23`)

## Function to calculate average statistics per gene number with Standard Error
calculate_average_statistics_per_gene_number <- function(data_list) {
  # Loop through each distance (name of the list element)
  result_list <- lapply(names(data_list), function(num_name) {
    df <- data_list[[num_name]]
    
    # Calculate sample sizes (excluding NAs)
    n_mean_ratio <- sum(!is.na(df$mean_ratio))
    n_p_value    <- sum(!is.na(df$p_value))
    
    # Calculate average and standard error for mean_ratio and p_value
    avg_mean_ratio <- mean(df$mean_ratio, na.rm = TRUE)
    se_mean_ratio  <- sd(df$mean_ratio, na.rm = TRUE) / sqrt(n_mean_ratio)
    
    avg_p_value <- mean(df$p_value, na.rm = TRUE)
    se_p_value  <- sd(df$p_value, na.rm = TRUE) / sqrt(n_p_value)
    
    # Create a row with all statistics
    data.frame(
      gene_number = as.numeric(num_name),
      average_mean_ratio = avg_mean_ratio,
      SE_mean_ratio = se_mean_ratio,
      average_p_value = avg_p_value,
      SE_p_value = se_p_value
    )
  })
  
  # Combine all rows into a single dataframe
  result_df <- do.call(rbind, result_list)
  rownames(result_df) <- NULL
  return(result_df)
}

##Apply the calculate_average_statistics_per_distance function for SRP066994, SRP048222, SRP048971
SRP066994_differing_num_genes_statistics_averages <- calculate_average_statistics_per_gene_number(SRP066994_differing_num_genes_means_ratios_and_pvalues_list)
SRP048222_differing_num_genes_statistics_averages <- calculate_average_statistics_per_gene_number(SRP048222_differing_num_genes_means_ratios_and_pvalues_list)
SRP048971_differing_num_genes_statistics_averages <- calculate_average_statistics_per_gene_number(SRP048971_differing_num_genes_means_ratios_and_pvalues_list)

##Checkpoint 20
#view(SRP066994_differing_num_genes_statistics_averages)
#view(SRP048222_differing_num_genes_statistics_averages)
#view(SRP048971_differing_num_genes_statistics_averages)

# Function to plot average mean ratios and p-values across gene numbers with SE bands
plot_average_mean_ratios_and_p_values_gene_numbers <- function(df) {
  # Plot 1: average_mean_ratio vs gene number with SE bands
  p1 <- ggplot(df, aes(x = gene_number, y = average_mean_ratio)) +
    geom_line(color = "blue") +
    geom_point(color = "blue") +
    geom_ribbon(aes(ymin = average_mean_ratio - SE_mean_ratio,
                    ymax = average_mean_ratio + SE_mean_ratio),
                alpha = 0.2, fill = "blue") +
    labs(title = "Average Mean Ratio vs Gene Number",
         x = "Gene Number",
         y = "Average Mean Ratio") +
    theme_minimal()
  
  # Plot 2: average_p_value vs gene number with SE bands
  p2 <- ggplot(df, aes(x = gene_number, y = average_p_value)) +
    geom_line(color = "red") +
    geom_point(color = "red") +
    geom_ribbon(aes(ymin = average_p_value - SE_p_value,
                    ymax = average_p_value + SE_p_value),
                alpha = 0.2, fill = "red") +
    labs(title = "Average p-Value vs Gene Number",
         x = "Gene Number",
         y = "Average p-Value") +
    theme_minimal()
  
  return(list(mean_ratio_plot = p1, p_value_plot = p2))
}

##Apply the plot_average_mean_ratios_and_p_values for SRP066994, SRP048222, SRP048971
SRP066994_gene_numbers_plots_list<- plot_average_mean_ratios_and_p_values_gene_numbers(SRP066994_differing_num_genes_statistics_averages)
SRP048222_gene_numbers_plots_list <- plot_average_mean_ratios_and_p_values_gene_numbers(SRP048222_differing_num_genes_statistics_averages)
SRP048971_gene_numbers_plots_list <- plot_average_mean_ratios_and_p_values_gene_numbers(SRP048971_differing_num_genes_statistics_averages)

##Checkpoint 21
#print(SRP066994_gene_numbers_plots_list$mean_ratio_plot)
#print(SRP066994_gene_numbers_plots_list$p_value_plot)
#print(SRP048222_gene_numbers_plots_list$mean_ratio_plot)
#print(SRP048222_gene_numbers_plots_list$p_value_plot)
print(SRP048971_gene_numbers_plots_list$mean_ratio_plot)
print(SRP048971_gene_numbers_plots_list$p_value_plot)
