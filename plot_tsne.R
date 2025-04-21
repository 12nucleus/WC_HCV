# --- Package Management ---
# List of required packages
required_packages <- c("dplyr", "tidyr", "stringr", "Rtsne", "ggplot2", "readr", "RColorBrewer") # Added RColorBrewer

# Check and install missing packages using a loop
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing missing package: ", pkg)
    tryCatch({
      install.packages(pkg, repos = "http://cran.us.r-project.org")
    }, error = function(e) {
      stop("Failed to install package '", pkg, "'. Please install manually.\nOriginal error: ", e$message)
    })
  }
}

# Load necessary libraries
message("Loading required libraries...")
library(dplyr)
library(tidyr)
library(stringr)
library(Rtsne)
library(ggplot2)
library(readr)
library(RColorBrewer) # Added RColorBrewer
message("Libraries loaded successfully.")
# --- End Package Management ---

# --- Configuration & Argument Parsing ---
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Usage: Rscript plot_tsne.R <input_pmatrix_file>", call. = FALSE)
} else {
  input_file <- args[1]
}

# Generate output filename based on input
output_plot <- sub("^(.*?)(\\.[^.]+)?$", "\\1_tsne.pdf", basename(input_file)) # Output as PDF

perplexity_value <- 15               # t-SNE perplexity (adjust based on data size, typically 5-50)
iterations <- 1000                   # t-SNE iterations
# --- End Configuration ---

# --- Data Loading and Preparation ---
message("Reading distance data from: ", input_file)
# Read the 3-column file, assuming no header and tab-separated
# Columns: Seq1, Seq2, Distance
tryCatch({
  # Use readr::read_tsv for more robust TSV parsing
  pairwise_distances <- readr::read_tsv(input_file, col_names = c("Seq1", "Seq2", "Distance"), col_types = "ccd", show_col_types = FALSE) # Changed last col_type to 'd' for double
}, error = function(e) {
  stop("Error reading input file '", input_file, "' using readr::read_tsv. Ensure it's a 3-column (Seq1, Seq2, Distance) tab-separated file with no header. Original error: ", e$message)
})

# Trim whitespace from sequence names
pairwise_distances <- pairwise_distances %>%
  mutate(
    Seq1 = trimws(Seq1),
    Seq2 = trimws(Seq2)
  )

# Check if any distances were read
if (nrow(pairwise_distances) == 0) {
    stop("No pairwise distances were read from the input file: ", input_file)
}

message("Extracting unique sequences...")
# Get all unique sequence names and trim whitespace
unique_sequences <- unique(trimws(c(pairwise_distances$Seq1, pairwise_distances$Seq2)))
n_seq <- length(unique_sequences)
message("Found ", n_seq, " unique sequences.")

if (n_seq < 2) {
    stop("Less than 2 unique sequences found. Cannot create distance matrix for t-SNE.")
}

message("Creating distance matrix...")
# Create an empty square matrix
dist_matrix <- matrix(NA, nrow = n_seq, ncol = n_seq, dimnames = list(unique_sequences, unique_sequences))

# Get the established dimension names for checking
matrix_dim_names <- rownames(dist_matrix) # or colnames, they are the same
missing_names_reported <- list() # To avoid flooding warnings

# Fill the matrix with distances
for (i in 1:nrow(pairwise_distances)) {
  # seq1 and seq2 already trimmed from mutate step above
  seq1 <- pairwise_distances$Seq1[i]
  seq2 <- pairwise_distances$Seq2[i]
  dist_val <- pairwise_distances$Distance[i]

  # Check if names exist in matrix dimensions before assigning
  seq1_exists <- seq1 %in% matrix_dim_names
  seq2_exists <- seq2 %in% matrix_dim_names

  if (seq1_exists && seq2_exists) {
      # Matrix is symmetric
      dist_matrix[seq1, seq2] <- dist_val
      dist_matrix[seq2, seq1] <- dist_val
  } else {
      if (!seq1_exists && !(seq1 %in% missing_names_reported)) {
          warning("Sequence name '", seq1, "' from input file not found in matrix dimensions. Skipping pair.")
          missing_names_reported[[seq1]] <- TRUE
      }
      if (!seq2_exists && !(seq2 %in% missing_names_reported)) {
          warning("Sequence name '", seq2, "' from input file not found in matrix dimensions. Skipping pair.")
          missing_names_reported[[seq2]] <- TRUE
      }
  }
}

# Set diagonal to 0 (distance to self)
diag(dist_matrix) <- 0

# Check for NAs (missing pairs) - Rtsne might handle this, or might need imputation
na_count <- sum(is.na(dist_matrix))
if (na_count > 0) {
  warning("Distance matrix contains ", na_count, " NA values (missing pairs). Attempting to proceed, but t-SNE might fail. Consider imputing missing values (e.g., with max distance) if needed.")
  # Optional: Impute NAs with max observed distance or a large value
  max_dist <- max(dist_matrix, na.rm = TRUE)
  if (is.finite(max_dist)) {
      dist_matrix[is.na(dist_matrix)] <- max_dist * 1.1 # Or some other strategy
      warning("Imputed NAs with ", max_dist * 1.1)
  } else {
      warning("Could not determine max distance to impute NAs. Proceeding with NAs.")
  }
}

# --- Metadata Extraction ---
message("Extracting metadata (Sample Name, Size)...")
metadata <- data.frame(SequenceName = unique_sequences, stringsAsFactors = FALSE) %>%
  mutate(
    # Extract Sample Name (look for GPxxxx or BL-xxxx pattern first)
    SampleName = str_extract(SequenceName, "GP\\d+|BL-[A-Za-z0-9]+"),
    # Fallback: If NA, try extracting content between first and last underscore (assuming number_name_number format)
    SampleName = ifelse(is.na(SampleName), str_extract(SequenceName, "(?<=^\\d+_)[^_]+(?=_\\d+$)"), SampleName),
    # Extract Size (the number at the very end)
    SequenceSize_str = str_extract(SequenceName, "\\d+$"),
    SequenceSize = suppressWarnings(as.numeric(SequenceSize_str)) # Suppress warnings for NAs introduced by coercion
  ) %>%
  select(SequenceName, SampleName, SequenceSize)

# Check if extraction worked
if(any(is.na(metadata$SampleName)) || any(is.na(metadata$SequenceSize))) {
  warning("Could not extract SampleName or SequenceSize for some sequences. Check the sequence naming format and the regex pattern.")
  # Print problematic names
  # print("Problematic sequence names:")
  # print(metadata[is.na(metadata$SampleName) | is.na(metadata$SequenceSize), ])
}
# Ensure SampleSize is numeric and handle potential NAs from conversion
metadata <- metadata %>% mutate(SequenceSize = ifelse(is.na(SequenceSize), 1, SequenceSize)) # Default size 1 if NA

# --- Run t-SNE ---
# Adjust perplexity if needed based on number of sequences
effective_perplexity <- perplexity_value
if (n_seq - 1 < 3 * perplexity_value) {
    effective_perplexity <- max(1, floor((n_seq - 1) / 3))
    warning("Perplexity is too high for the number of samples (", n_seq, "). Adjusting perplexity to ", effective_perplexity)
}
effective_perplexity <- min(effective_perplexity, n_seq - 1) # Perplexity must be less than N
effective_perplexity <- max(1, effective_perplexity) # Ensure perplexity is at least 1

message("Running t-SNE (perplexity=", effective_perplexity, ", iterations=", iterations, ")...")
# Set seed for reproducibility
set.seed(42)
# Run Rtsne directly on the matrix, indicating it's a distance matrix
# Use check_duplicates=FALSE if sequence names might appear identical but represent distinct points
tsne_results <- tryCatch({
    Rtsne(dist_matrix, is_distance = TRUE, perplexity = effective_perplexity, max_iter = iterations, check_duplicates = FALSE, pca = FALSE) # No PCA needed for distance matrix
}, error = function(e) {
    warning("Rtsne failed: ", e$message)
    return(NULL)
})

# --- Prepare Data for Plotting ---
if (!is.null(tsne_results)) {
    message("Preparing data for plotting...")
    tsne_data <- data.frame(
      SequenceName = unique_sequences, # Ensure order matches dist_matrix rows
      TSNE1 = tsne_results$Y[, 1],
      TSNE2 = tsne_results$Y[, 2]
    )

    # Merge with metadata
    plot_data <- left_join(tsne_data, metadata, by = "SequenceName")

    # Ensure SampleName is treated as a factor for coloring and legend
    plot_data <- plot_data %>%
      mutate(
        SampleName = factor(SampleName), # Convert SampleName to factor
        # Calculate log size, handle size 0 or 1 appropriately (log(1)=0)
        LogSize = log2(SequenceSize + 1) # Add 1 to avoid log(0) issues
      )

    # --- Debugging: Inspect plot_data ---
    message("Structure of plot_data before plotting:")
    print(str(plot_data))
    message("First few rows of plot_data:")
    print(head(plot_data))
    message("Unique Sample Names (Factor Levels) being used:")
    print(levels(plot_data$SampleName))
    # --- End Debugging ---

    # Generate a color palette for all unique samples
    unique_sample_names <- levels(plot_data$SampleName)
    num_samples <- length(unique_sample_names)
    message("Found ", num_samples, " unique samples for coloring: ", paste(unique_sample_names, collapse=", "))

    # Use a Brewer palette, recycling if necessary
    if (num_samples <= 8) {
        # Use a qualitative palette suitable for fewer categories
        palette_name <- "Set1"
        max_colors <- brewer.pal.info[palette_name, "maxcolors"]
        if (num_samples > max_colors) { palette_name <- "Set2"; max_colors <- brewer.pal.info[palette_name, "maxcolors"] }
        if (num_samples > max_colors) { palette_name <- "Paired"; max_colors <- brewer.pal.info[palette_name, "maxcolors"] }
        if (num_samples > max_colors) { palette_name <- "Dark2"; max_colors <- brewer.pal.info[palette_name, "maxcolors"] }

        # Fallback if still too many
        if (num_samples > max_colors) {
             colors_to_use <- colorRampPalette(brewer.pal(max_colors, palette_name))(num_samples)
        } else {
             colors_to_use <- brewer.pal(num_samples, palette_name)
        }

    } else {
        # Use a palette generator for more colors if needed
        # Use a palette known for many distinct colors
        qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
        col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
        if (num_samples > length(col_vector)) {
             # If still not enough, recycle with slight modification or use rainbow
             colors_to_use <- colorRampPalette(col_vector)(num_samples)
             # colors_to_use <- rainbow(num_samples)
        } else {
             colors_to_use <- col_vector[1:num_samples]
        }
    }

    # Create the named color map
    color_map <- setNames(colors_to_use, unique_sample_names)

    # --- Prepare Edge Data (Links <= 10 SNPs) ---
    message("Preparing edge data for plotting (Distance <= 10)...")
    edge_data_raw <- pairwise_distances %>%
      filter(Distance <= 10 & Seq1 != Seq2) # Filter links and remove self-loops

    # Ensure only one edge per pair of sequences, regardless of direction
    edge_data_unique <- edge_data_raw %>%
      rowwise() %>%
      mutate(pair_key = paste(sort(c(Seq1, Seq2)), collapse = "___")) %>% # Create unique key for the pair
      ungroup() %>%
      distinct(pair_key, .keep_all = TRUE) %>% # Keep only one row per unique pair
      select(-pair_key) # Remove the temporary key

    # Join to get coordinates and sample names for edges
    edge_data_coords <- edge_data_unique %>%
      # Join for Seq1
      left_join(plot_data %>% select(SequenceName, SampleName, TSNE1, TSNE2), by = c("Seq1" = "SequenceName")) %>%
      rename(SampleName1 = SampleName, x = TSNE1, y = TSNE2) %>%
      # Join for Seq2
      left_join(plot_data %>% select(SequenceName, SampleName, TSNE1, TSNE2), by = c("Seq2" = "SequenceName")) %>%
      rename(SampleName2 = SampleName, xend = TSNE1, yend = TSNE2) %>%
      # Filter out edges within the same sample and edges missing coordinates
      filter(SampleName1 != SampleName2 & !is.na(x) & !is.na(y) & !is.na(xend) & !is.na(yend))

    # Final edge data for plotting
    edge_data <- edge_data_coords

    message("Found ", nrow(edge_data), " unique inter-sample edges with distance <= 10 to plot.")

    # --- Calculate Dynamic Alpha for Edges ---
    n_edges <- nrow(edge_data)
    max_alpha <- 0.3
    min_alpha <- 0.005
    target_edges_for_max_alpha <- 25 # Number of edges where alpha is max_alpha
    
    # Calculate alpha based on inverse relationship, clamped between min/max
    if (n_edges > 0) {
      dynamic_alpha <- (target_edges_for_max_alpha / n_edges) * max_alpha
      dynamic_alpha <- max(min_alpha, min(max_alpha, dynamic_alpha)) # Clamp between min and max
    } else {
      dynamic_alpha <- max_alpha # Default if no edges
    }
    message("Using dynamic alpha for edges: ", round(dynamic_alpha, 3))

    # --- Generate Plot ---
    message("Generating plot...")
    p <- ggplot(plot_data, aes(x = TSNE1, y = TSNE2, color = SampleName, size = LogSize)) + # Use SampleName directly
      # Add segments first so points are drawn on top
      geom_segment(data = edge_data, aes(x=x, y=y, xend=xend, yend=yend),
                   color = "black", alpha = dynamic_alpha, linewidth = 0.1, inherit.aes = FALSE) + # Use dynamic alpha
      geom_point(alpha = 0.6) + # Slightly increase alpha for points
      scale_color_manual(
          name = "Sample", # Legend title for color
          values = color_map, # Use the generated color map for all samples
          labels = unique_sample_names # Ensure labels match the factor levels
          ) +
      scale_size_continuous(name = "Log2") + # Legend title for size
      labs(
        #title = "t-SNE Visualization of Sequence Distances", # Title can be added if desired
        subtitle = paste("Edges connect sequences from different samples with SNP distance <= 10"), # Updated subtitle
        x = "t-SNE Dimension 1",
        y = "t-SNE Dimension 2"
      ) +
      theme_light() + # Theme with light grey background and white grid
      theme(legend.position = "right")

    # --- Save Plot ---
    message("Saving plot to: ", output_plot)
    ggsave(output_plot, plot = p, device = "pdf", width = 10, height = 8) # Save as PDF

    message("Script finished successfully.")
} else {
    message("Skipping plot generation because Rtsne failed.")
}