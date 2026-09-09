# ============================================================
# dyngen strong cluster simulation
#
# 5 classes:
#   genes = 600
#   markers/class = 100
#   cells = 200, 400, 600, 800, 1000
#
# 10 classes:
#   genes = 1100
#   markers/class = 100
#   cells = 200, 200, 300, 300, 500,
#           500, 800, 1000, 1000, 1200
#
# Marker strength:
#   transcription_rate fold ~ Uniform(8, 12)
#
# Output:
#   spliced.mtx       genes x cells
#   unspliced.mtx     genes x cells
#   subclass_label
#   barcode
#   gene_name
#
# Maximum:
#   max(S + U) <= 800
# ============================================================


suppressPackageStartupMessages({
  library(dyngen)
  library(Matrix)
})


# ============================================================
# 1. Global settings
# ============================================================

OUT_DIR <- "dyngen_simulation"

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


CACHE_DIR <- tools::R_user_dir(
  "dyngen",
  "data"
)

dir.create(
  CACHE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


N_CORES <- parallel::detectCores(
  logical = TRUE
)

if (is.na(N_CORES)) {
  N_CORES <- 1L
}

N_CORES <- max(
  1L,
  N_CORES - 1L
)

options(Ncpus = N_CORES)

options(
  dyngen_download_cache_dir =
    CACHE_DIR
)



# ============================================================
# 2. Main parameters
# ============================================================

MARKERS_PER_CLASS <- 100L


# ------------------------------------------------------------
# Marker strength
# ------------------------------------------------------------

MARKER_FOLD_MIN <- 8.0

MARKER_FOLD_MAX <- 12.0


# ------------------------------------------------------------
# Maximum single gene count
#
# S + U <= 800
# ------------------------------------------------------------

MAX_COUNT <- 800L


# ------------------------------------------------------------
# Small kinetic variation
# ------------------------------------------------------------

KINETIC_NOISE_SD <- 0.001


# ------------------------------------------------------------
# SSA
# ------------------------------------------------------------

SSA_TAU <- 0.02

CENSUS_INTERVAL <- 2

SIM_TIME_MULTIPLIER <- 2


# More simulated census points than final cells
CENSUS_RATIO <- 6


# ------------------------------------------------------------
# Library-size variation
# ------------------------------------------------------------

LIB_CV_LOG <- 0.12


# Medium-high expression intensity
LIB_MEDIAN_5CLASS <- 30000

LIB_MEDIAN_10CLASS <- 52000


# ------------------------------------------------------------
# Seeds
# ------------------------------------------------------------

SEED_5CLASS <- 20260805L

SEED_10CLASS <- 20260810L



# ============================================================
# 3. Construct library-size reference
# ============================================================

make_library_reference <- function(
    median_library,
    n_reference = 5000L,
    seed = 1L
) {
  
  set.seed(seed)
  
  libs <- round(
    exp(
      rnorm(
        n_reference,
        mean = log(median_library),
        sd = LIB_CV_LOG
      )
    )
  )
  
  
  lower <- round(
    median_library * 0.70
  )
  
  upper <- round(
    median_library * 1.30
  )
  
  
  libs <- pmax(
    lower,
    pmin(
      upper,
      libs
    )
  )
  
  
  Matrix::sparseMatrix(
    
    i = seq_len(
      n_reference
    ),
    
    j = rep(
      1L,
      n_reference
    ),
    
    x = libs,
    
    dims = c(
      n_reference,
      1L
    )
  )
}



# ============================================================
# 4. Determine number of independent simulations
# ============================================================

calculate_num_simulations <- function(
    n_cells,
    total_time,
    census_interval
) {
  
  census_per_sim <- max(
    1L,
    floor(
      total_time /
        census_interval
    )
  )
  
  
  n_sim <- ceiling(
    n_cells *
      CENSUS_RATIO /
      census_per_sim
  )
  
  
  max(
    10L,
    as.integer(
      n_sim
    )
  )
}



# ============================================================
# 5. Enforce max(S+U) <= 800
#
# Whole-cell binomial thinning is used rather than simply
# truncating individual count values.
# ============================================================

thin_cell_to_cap <- function(
    U,
    S,
    cap = 800L,
    seed = NULL
) {
  
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  
  stopifnot(
    identical(
      dim(U),
      dim(S)
    )
  )
  
  
  for (i in seq_len(nrow(U))) {
    
    repeat {
      
      total_i <- U[i, ] + S[i, ]
      
      current_max <- max(
        total_i
      )
      
      
      if (current_max <= cap) {
        break
      }
      
      
      # Leave a little margin below 800
      target <- cap * 0.94
      
      
      p <- target /
        current_max
      
      
      p <- min(
        0.995,
        p
      )
      
      p <- max(
        0.01,
        p
      )
      
      
      U[i, ] <- rbinom(
        
        n = ncol(U),
        
        size = as.integer(
          round(
            U[i, ]
          )
        ),
        
        prob = p
      )
      
      
      S[i, ] <- rbinom(
        
        n = ncol(S),
        
        size = as.integer(
          round(
            S[i, ]
          )
        ),
        
        prob = p
      )
    }
  }
  
  
  list(
    U = U,
    S = S
  )
}



# ============================================================
# 6. Main simulation function
# ============================================================

simulate_cluster_dataset <- function(
    n_classes,
    n_genes,
    class_cell_numbers,
    output_prefix,
    median_library,
    seed
) {
  
  cat("\n")
  cat("============================================================\n")
  cat("Dataset:", output_prefix, "\n")
  cat("Classes:", n_classes, "\n")
  cat("Genes:", n_genes, "\n")
  cat("Cells:", sum(class_cell_numbers), "\n")
  cat("Markers/class:", MARKERS_PER_CLASS, "\n")
  cat(
    "Marker kinetic fold:",
    MARKER_FOLD_MIN,
    "-",
    MARKER_FOLD_MAX,
    "\n"
  )
  cat("============================================================\n\n")
  
  
  stopifnot(
    length(class_cell_numbers) ==
      n_classes
  )
  
  
  n_marker_total <-
    n_classes *
    MARKERS_PER_CLASS
  
  
  stopifnot(
    n_marker_total <
      n_genes
  )
  
  
  # ==========================================================
  # 6.1 Common backbone
  # ==========================================================
  
  set.seed(seed)
  
  
  backbone <-
    backbone_linear_simple()
  
  
  NUM_TFS <-
    nrow(
      backbone$module_info
    )
  
  
  # ==========================================================
  # Marker genes
  #
  # 5-class:
  #   500 marker genes
  #   100 shared genes
  #
  # 10-class:
  #   1000 marker genes
  #   100 shared genes
  # ==========================================================
  
  NUM_HKS <-
    n_marker_total
  
  
  NUM_TARGETS <-
    n_genes -
    NUM_HKS -
    NUM_TFS
  
  
  if (NUM_TARGETS < 1L) {
    
    stop(
      "Not enough genes for common GRN."
    )
  }
  
  
  cat(
    "TF genes:",
    NUM_TFS,
    "\n"
  )
  
  cat(
    "Marker-carrier genes:",
    NUM_HKS,
    "\n"
  )
  
  cat(
    "Shared target genes:",
    NUM_TARGETS,
    "\n"
  )
  
  
  # ==========================================================
  # Simulation time
  # ==========================================================
  
  total_time <-
    simtime_from_backbone(
      backbone
    ) *
    SIM_TIME_MULTIPLIER
  
  
  max_cells <-
    max(
      class_cell_numbers
    )
  
  
  n_sim_initial <-
    calculate_num_simulations(
      
      n_cells =
        max_cells,
      
      total_time =
        total_time,
      
      census_interval =
        CENSUS_INTERVAL
    )
  
  
  # ==========================================================
  # 6.2 Desired library-size distribution
  # ==========================================================
  
  reference_counts <-
    make_library_reference(
      
      median_library =
        median_library,
      
      n_reference =
        5000L,
      
      seed =
        seed + 10L
    )
  
  
  # ==========================================================
  # 6.3 Initialise dyngen
  # ==========================================================
  
  config <-
    initialise_model(
      
      backbone =
        backbone,
      
      num_cells =
        max_cells,
      
      num_tfs =
        NUM_TFS,
      
      num_targets =
        NUM_TARGETS,
      
      num_hks =
        NUM_HKS,
      
      distance_metric =
        "euclidean",
      
      
      kinetics_params =
        kinetics_default(),
      
      
      simulation_params =
        simulation_default(
          
          total_time =
            total_time,
          
          census_interval =
            CENSUS_INTERVAL,
          
          ssa_algorithm =
            ssa_etl(
              tau = SSA_TAU
            ),
          
          experiment_params =
            simulation_type_wild_type(
              num_simulations =
                n_sim_initial
            ),
          
          compute_dimred =
            FALSE,
          
          compute_rna_velocity =
            FALSE,
          
          kinetics_noise_function =
            kinetics_noise_simple(
              mean = 1,
              sd = KINETIC_NOISE_SD
            )
        ),
      
      
      experiment_params =
        experiment_snapshot(
          
          realcount =
            reference_counts,
          
          map_reference_cpm =
            FALSE,
          
          map_reference_ls =
            TRUE,
          
          weight_bw =
            0.1
        ),
      
      
      verbose =
        TRUE,
      
      download_cache_dir =
        CACHE_DIR,
      
      num_cores =
        N_CORES
    )
  
  
  # ==========================================================
  # 6.4 Generate common GRN and kinetics
  # ==========================================================
  
  set.seed(
    seed + 20L
  )
  
  
  base_model <-
    config |>
    generate_tf_network() |>
    generate_feature_network() |>
    generate_kinetics() |>
    generate_gold_standard()
  
  
  if (
    nrow(
      base_model$feature_info
    ) != n_genes
  ) {
    
    stop(
      "Expected ",
      n_genes,
      " genes, but dyngen returned ",
      nrow(
        base_model$feature_info
      )
    )
  }
  
  
  # ==========================================================
  # 6.5 Find marker-carrier genes
  # ==========================================================
  
  hk_ids <-
    base_model$feature_info$feature_id[
      base_model$feature_info$is_hk %in% TRUE
    ]
  
  
  if (
    length(hk_ids) !=
    n_marker_total
  ) {
    
    stop(
      "Expected ",
      n_marker_total,
      " marker-carrier genes, found ",
      length(hk_ids)
    )
  }
  
  
  # Random ordering but reproducible
  set.seed(
    seed + 30L
  )
  
  
  marker_ids_flat <-
    sample(
      hk_ids,
      size =
        n_marker_total,
      replace =
        FALSE
    )
  
  
  marker_sets <-
    split(
      
      marker_ids_flat,
      
      rep(
        seq_len(
          n_classes
        ),
        each =
          MARKERS_PER_CLASS
      )
    )
  
  
  # ==========================================================
  # 6.6 Final gene ordering
  #
  # Class1:
  # Gene1-Gene100
  #
  # Class2:
  # Gene101-Gene200
  #
  # ...
  # ==========================================================
  
  all_feature_ids <-
    base_model$feature_info$feature_id
  
  
  shared_ids <-
    setdiff(
      all_feature_ids,
      marker_ids_flat
    )
  
  
  gene_order_ids <-
    c(
      marker_ids_flat,
      shared_ids
    )
  
  
  stopifnot(
    length(
      gene_order_ids
    ) ==
      n_genes
  )
  
  
  gene_names <-
    paste0(
      "Gene",
      seq_len(
        n_genes
      )
    )
  
  
  marker_class <-
    c(
      
      rep(
        paste0(
          "Class",
          seq_len(
            n_classes
          )
        ),
        each =
          MARKERS_PER_CLASS
      ),
      
      rep(
        "shared",
        n_genes -
          n_marker_total
      )
    )
  
  
  # ==========================================================
  # 6.7 Each marker gets a different 8-12 fold value
  # ==========================================================
  
  set.seed(
    seed + 40L
  )
  
  
  marker_fold_values <-
    runif(
      
      n =
        n_marker_total,
      
      min =
        MARKER_FOLD_MIN,
      
      max =
        MARKER_FOLD_MAX
    )
  
  
  names(
    marker_fold_values
  ) <-
    marker_ids_flat
  
  
  cat(
    "\nActual marker kinetic-fold range:",
    round(
      min(marker_fold_values),
      3
    ),
    "-",
    round(
      max(marker_fold_values),
      3
    ),
    "\n"
  )
  
  
  cat(
    "Mean marker kinetic fold:",
    round(
      mean(marker_fold_values),
      3
    ),
    "\n\n"
  )
  
  
  # ==========================================================
  # Containers
  # ==========================================================
  
  S_list <-
    vector(
      "list",
      n_classes
    )
  
  
  U_list <-
    vector(
      "list",
      n_classes
    )
  
  
  label_list <-
    vector(
      "list",
      n_classes
    )
  
  
  barcode_list <-
    vector(
      "list",
      n_classes
    )
  
  
  # ==========================================================
  # 6.8 Generate each class
  # ==========================================================
  
  for (k in seq_len(n_classes)) {
    
    n_cells_k <-
      class_cell_numbers[k]
    
    
    cat("\n")
    cat("------------------------------------------------------------\n")
    cat("Class:", k, "\n")
    cat("Cells:", n_cells_k, "\n")
    cat("Markers:", MARKERS_PER_CLASS, "\n")
    cat("------------------------------------------------------------\n")
    
    
    # Same common model
    model_k <-
      base_model
    
    
    model_k$numbers$num_cells <-
      as.integer(
        n_cells_k
      )
    
    
    # ========================================================
    # Number of SSA simulations
    # ========================================================
    
    n_sim_k <-
      calculate_num_simulations(
        
        n_cells =
          n_cells_k,
        
        total_time =
          total_time,
        
        census_interval =
          CENSUS_INTERVAL
      )
    
    
    cat(
      "SSA simulations:",
      n_sim_k,
      "\n"
    )
    
    
    model_k$simulation_params$experiment_params <-
      simulation_type_wild_type(
        num_simulations =
          n_sim_k
      )
    
    
    model_k$experiment_params <-
      experiment_snapshot(
        
        realcount =
          reference_counts,
        
        map_reference_cpm =
          FALSE,
        
        map_reference_ls =
          TRUE,
        
        weight_bw =
          0.1
      )
    
    
    # ========================================================
    # Class-specific marker enhancement
    #
    # transcription_rate(class)
    #
    # =
    #
    # baseline transcription_rate
    #
    # ×
    #
    # Uniform(8,12)
    # ========================================================
    
    markers_k <-
      marker_sets[[k]]
    
    
    marker_idx <-
      match(
        
        markers_k,
        
        model_k$feature_info$feature_id
      )
    
    
    if (
      any(
        is.na(
          marker_idx
        )
      )
    ) {
      
      stop(
        "Marker lookup failed for Class ",
        k
      )
    }
    
    
    folds_k <-
      marker_fold_values[
        markers_k
      ]
    
    
    original_rate <-
      model_k$feature_info$transcription_rate[
        marker_idx
      ]
    
    
    model_k$feature_info$transcription_rate[
      marker_idx
    ] <-
      original_rate *
      folds_k
    
    
    cat(
      "Marker fold mean:",
      round(
        mean(
          folds_k
        ),
        3
      ),
      "\n"
    )
    
    
    cat(
      "Marker fold range:",
      round(
        min(
          folds_k
        ),
        3
      ),
      "-",
      round(
        max(
          folds_k
        ),
        3
      ),
      "\n"
    )
    
    
    # ========================================================
    # Generate cells and transcript sampling
    # ========================================================
    
    set.seed(
      seed +
        1000L +
        k
    )
    
    
    model_k <-
      model_k |>
      generate_cells() |>
      generate_experiment()
    
    
    # ========================================================
    # Extract:
    #
    # U = pre-mRNA
    # S = mature mRNA
    #
    # Current:
    # cells x genes
    # ========================================================
    
    U_k <-
      as.matrix(
        
        model_k$experiment$counts_premrna[
          ,
          gene_order_ids,
          drop = FALSE
        ]
      )
    
    
    S_k <-
      as.matrix(
        
        model_k$experiment$counts_mrna[
          ,
          gene_order_ids,
          drop = FALSE
        ]
      )
    
    
    if (
      nrow(S_k) !=
      n_cells_k
    ) {
      
      stop(
        "Class ",
        k,
        ": expected ",
        n_cells_k,
        " cells, got ",
        nrow(S_k)
      )
    }
    
    
    # ========================================================
    # Strict max count <= 800
    # ========================================================
    
    cap_result <-
      thin_cell_to_cap(
        
        U =
          U_k,
        
        S =
          S_k,
        
        cap =
          MAX_COUNT,
        
        seed =
          seed +
          2000L +
          k
      )
    
    
    U_k <-
      cap_result$U
    
    
    S_k <-
      cap_result$S
    
    
    stopifnot(
      max(
        U_k +
          S_k
      ) <=
        MAX_COUNT
    )
    
    
    # ========================================================
    # Labels
    # ========================================================
    
    label_k <-
      rep(
        paste0(
          "Class",
          k
        ),
        n_cells_k
      )
    
    
    barcode_k <-
      sprintf(
        
        "%s_Class%02d_Cell%05d",
        
        output_prefix,
        
        k,
        
        seq_len(
          n_cells_k
        )
      )
    
    
    S_list[[k]] <-
      S_k
    
    
    U_list[[k]] <-
      U_k
    
    
    label_list[[k]] <-
      label_k
    
    
    barcode_list[[k]] <-
      barcode_k
    
    
    cat(
      "Median library size:",
      round(
        median(
          rowSums(
            S_k +
              U_k
          )
        )
      ),
      "\n"
    )
    
    
    cat(
      "Maximum total count:",
      max(
        S_k +
          U_k
      ),
      "\n"
    )
    
    
    rm(
      model_k,
      S_k,
      U_k,
      cap_result
    )
    
    gc()
  }
  
  
  # ==========================================================
  # 6.9 Merge classes
  #
  # Current:
  # cells x genes
  # ==========================================================
  
  S_all <-
    do.call(
      rbind,
      S_list
    )
  
  
  U_all <-
    do.call(
      rbind,
      U_list
    )
  
  
  labels <-
    unlist(
      label_list,
      use.names = FALSE
    )
  
  
  barcodes <-
    unlist(
      barcode_list,
      use.names = FALSE
    )
  
  
  stopifnot(
    
    nrow(S_all) ==
      sum(
        class_cell_numbers
      ),
    
    ncol(S_all) ==
      n_genes,
    
    identical(
      dim(S_all),
      dim(U_all)
    )
  )
  
  
  colnames(S_all) <-
    gene_names
  
  
  colnames(U_all) <-
    gene_names
  
  
  rownames(S_all) <-
    barcodes
  
  
  rownames(U_all) <-
    barcodes
  
  
  X_all <-
    S_all +
    U_all
  
  
  stopifnot(
    max(
      X_all
    ) <=
      MAX_COUNT
  )
  
  
  # ==========================================================
  # 6.10 Dataset summary
  # ==========================================================
  
  cat("\n")
  cat("============================================================\n")
  cat("FINAL DATA SUMMARY\n")
  cat("============================================================\n")
  
  
  cat(
    "Dataset:",
    output_prefix,
    "\n"
  )
  
  
  cat(
    "Shape cells x genes:",
    nrow(X_all),
    "x",
    ncol(X_all),
    "\n"
  )
  
  
  cat(
    "Maximum total count:",
    max(X_all),
    "\n"
  )
  
  
  cat(
    "Median library size:",
    round(
      median(
        rowSums(
          X_all
        )
      )
    ),
    "\n"
  )
  
  
  cat(
    "Mean library size:",
    round(
      mean(
        rowSums(
          X_all
        )
      )
    ),
    "\n"
  )
  
  
  # ==========================================================
  # 6.11 Marker separation diagnostic
  # ==========================================================
  
  marker_summary <-
    vector(
      "list",
      n_classes
    )
  
  
  for (
    k in seq_len(
      n_classes
    )
  ) {
    
    cell_idx <-
      labels ==
      paste0(
        "Class",
        k
      )
    
    
    marker_idx <-
      (
        (k - 1L) *
          MARKERS_PER_CLASS +
          1L
      ):(
        k *
          MARKERS_PER_CLASS
      )
    
    
    marker_mean_own <-
      mean(
        X_all[
          cell_idx,
          marker_idx,
          drop = FALSE
        ]
      )
    
    
    marker_mean_other <-
      mean(
        X_all[
          !cell_idx,
          marker_idx,
          drop = FALSE
        ]
      )
    
    
    marker_ratio <-
      marker_mean_own /
      (
        marker_mean_other +
          1e-8
      )
    
    
    folds_this_class <-
      marker_fold_values[
        marker_sets[[k]]
      ]
    
    
    marker_summary[[k]] <-
      data.frame(
        
        class =
          paste0(
            "Class",
            k
          ),
        
        n_cells =
          sum(
            cell_idx
          ),
        
        marker_start =
          min(
            marker_idx
          ),
        
        marker_end =
          max(
            marker_idx
          ),
        
        kinetic_fold_mean =
          mean(
            folds_this_class
          ),
        
        kinetic_fold_min =
          min(
            folds_this_class
          ),
        
        kinetic_fold_max =
          max(
            folds_this_class
          ),
        
        marker_mean_own =
          marker_mean_own,
        
        marker_mean_other =
          marker_mean_other,
        
        own_vs_other_ratio =
          marker_ratio
      )
    
    
    cat(
      "Class",
      k,
      "marker own/other ratio =",
      round(
        marker_ratio,
        3
      ),
      "\n"
    )
  }
  
  
  marker_summary <-
    do.call(
      rbind,
      marker_summary
    )
  
  
  # ==========================================================
  # 6.12 Marker information
  # ==========================================================
  
  marker_fold_out <-
    c(
      
      as.numeric(
        marker_fold_values[
          marker_ids_flat
        ]
      ),
      
      rep(
        1,
        n_genes -
          n_marker_total
      )
    )
  
  
  marker_table <-
    data.frame(
      
      gene_name =
        gene_names,
      
      marker_class =
        marker_class,
      
      kinetic_fold =
        marker_fold_out
    )
  
  
  write.table(
    
    marker_table,
    
    file =
      file.path(
        OUT_DIR,
        paste0(
          output_prefix,
          "_marker_table.tsv"
        )
      ),
    
    sep = "\t",
    
    quote = FALSE,
    
    row.names = FALSE
  )
  
  
  write.table(
    
    marker_summary,
    
    file =
      file.path(
        OUT_DIR,
        paste0(
          output_prefix,
          "_marker_summary.tsv"
        )
      ),
    
    sep = "\t",
    
    quote = FALSE,
    
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 6.13 Export genes x cells
  # ==========================================================
  
  S_sparse <-
    Matrix::Matrix(
      t(S_all),
      sparse = TRUE
    )
  
  
  U_sparse <-
    Matrix::Matrix(
      t(U_all),
      sparse = TRUE
    )
  
  
  # ==========================================================
  # spliced
  # ==========================================================
  
  Matrix::writeMM(
    
    S_sparse,
    
    file.path(
      OUT_DIR,
      paste0(
        output_prefix,
        "_spliced.mtx"
      )
    )
  )
  
  
  # ==========================================================
  # unspliced
  # ==========================================================
  
  Matrix::writeMM(
    
    U_sparse,
    
    file.path(
      OUT_DIR,
      paste0(
        output_prefix,
        "_unspliced.mtx"
      )
    )
  )
  
  
  # ==========================================================
  # gene names
  # ==========================================================
  
  writeLines(
    
    gene_names,
    
    file.path(
      OUT_DIR,
      paste0(
        output_prefix,
        "_gene_name.txt"
      )
    )
  )
  
  
  # ==========================================================
  # labels
  # ==========================================================
  
  writeLines(
    
    labels,
    
    file.path(
      OUT_DIR,
      paste0(
        output_prefix,
        "_subclass_label.txt"
      )
    )
  )
  
  
  # ==========================================================
  # barcode
  # ==========================================================
  
  writeLines(
    
    barcodes,
    
    file.path(
      OUT_DIR,
      paste0(
        output_prefix,
        "_barcode.txt"
      )
    )
  )
  
  
  # ==========================================================
  # marker class
  # ==========================================================
  
  writeLines(
    
    marker_class,
    
    file.path(
      OUT_DIR,
      paste0(
        output_prefix,
        "_marker_class.txt"
      )
    )
  )
  
  
  # ==========================================================
  # Summary file
  # ==========================================================
  
  dataset_summary <-
    data.frame(
      
      dataset =
        output_prefix,
      
      n_classes =
        n_classes,
      
      n_genes =
        n_genes,
      
      n_cells =
        nrow(
          X_all
        ),
      
      markers_per_class =
        MARKERS_PER_CLASS,
      
      marker_fold_min =
        MARKER_FOLD_MIN,
      
      marker_fold_max =
        MARKER_FOLD_MAX,
      
      max_count =
        max(
          X_all
        ),
      
      median_library_size =
        median(
          rowSums(
            X_all
          )
        ),
      
      mean_library_size =
        mean(
          rowSums(
            X_all
          )
        )
    )
  
  
  write.table(
    
    dataset_summary,
    
    file =
      file.path(
        OUT_DIR,
        paste0(
          output_prefix,
          "_summary.tsv"
        )
      ),
    
    sep = "\t",
    
    quote = FALSE,
    
    row.names = FALSE
  )
  
  
  cat("\n")
  cat(
    "Saved:",
    output_prefix,
    "\n"
  )
  
  cat("============================================================\n")
  
  
  invisible(
    list(
      marker_summary =
        marker_summary,
      dataset_summary =
        dataset_summary
    )
  )
}



# ============================================================
# 7. Generate 5-class dataset
#
# Gene1-Gene100:
#   Class1 marker
#
# Gene101-Gene200:
#   Class2 marker
#
# Gene201-Gene300:
#   Class3 marker
#
# Gene301-Gene400:
#   Class4 marker
#
# Gene401-Gene500:
#   Class5 marker
#
# Gene501-Gene600:
#   shared
# ============================================================

result_5 <-
  simulate_cluster_dataset(
    
    n_classes =
      5L,
    
    n_genes =
      600L,
    
    class_cell_numbers =
      c(
        200L,
        400L,
        600L,
        800L,
        1000L
      ),
    
    output_prefix =
      "dyngen_5class_600genes",
    
    median_library =
      LIB_MEDIAN_5CLASS,
    
    seed =
      SEED_5CLASS
  )



# ============================================================
# 8. Generate 10-class dataset
#
# Gene1-Gene100:
#   Class1
#
# ...
#
# Gene901-Gene1000:
#   Class10
#
# Gene1001-Gene1100:
#   shared
# ============================================================

result_10 <-
  simulate_cluster_dataset(
    
    n_classes =
      10L,
    
    n_genes =
      1100L,
    
    class_cell_numbers =
      c(
        200L,
        200L,
        300L,
        300L,
        500L,
        500L,
        800L,
        1000L,
        1000L,
        1200L
      ),
    
    output_prefix =
      "dyngen_10class_1100genes",
    
    median_library =
      LIB_MEDIAN_10CLASS,
    
    seed =
      SEED_10CLASS
  )


cat("\n")
cat("============================================================\n")
cat("ALL DATASETS FINISHED\n")
cat("Marker kinetic strength = 8-12\n")
cat("Kinetic noise SD = 0.001\n")
cat("Maximum S+U <= 800\n")
cat("============================================================\n")