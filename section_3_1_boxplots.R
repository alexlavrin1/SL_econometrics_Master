# ==============================================================================
# Section 3.1 - NRMSE boxplots
#
# This script does not run the simulation. It reads the saved simulation output
# from sim_results.rds and creates the six-panel NRMSE boxplot figure.
# Keep this script and sim_results.rds in the same project folder.
# ==============================================================================

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required. Install it once with: ",
    "install.packages('ggplot2')"
  )
}

get_script_directory <- function() {
  file_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)

  if (length(file_argument) > 0) {
    script_path <- sub("^--file=", "", file_argument[1])
    return(dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE)))
  }

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    if (!is.null(frames[[i]]$ofile)) {
      return(dirname(normalizePath(
        frames[[i]]$ofile,
        winslash = "/",
        mustWork = TRUE
      )))
    }
  }

  getwd()
}

project_directory <- get_script_directory()
results_file <- file.path(project_directory, "sim_results.rds")

if (!file.exists(results_file)) {
  stop(
    "Cannot find 'sim_results.rds'. Put it in the project folder: ",
    project_directory
  )
}

sim <- readRDS(results_file)

if (!is.list(sim) || is.null(sim$results)) {
  stop("'sim_results.rds' does not contain the expected 'results' data frame.")
}

results <- sim$results
required_columns <- c("nrmse", "scenario", "mechanism", "method", "rep")
missing_columns <- setdiff(required_columns, names(results))

if (length(missing_columns) > 0) {
  stop(
    "The results data frame is missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Obtain one NRMSE value per replication and method by averaging over X1-X3.
per_rep <- aggregate(
  nrmse ~ scenario + mechanism + method + rep,
  data = results,
  FUN = mean
)

per_rep$mechanism <- factor(
  per_rep$mechanism,
  levels = c("MCAR", "MAR", "MNAR")
)

per_rep$scenario <- factor(
  per_rep$scenario,
  levels = c("A", "B"),
  labels = c("Scenario A", "Scenario B")
)

per_rep$method <- factor(
  per_rep$method,
  levels = c("hybrid", "unweighted", "mean"),
  labels = c("Hybrid kNN", "Standard kNN", "Mean")
)

p_nrmse <- ggplot2::ggplot(
  per_rep,
  ggplot2::aes(x = method, y = nrmse, fill = method)
) +
  ggplot2::geom_boxplot(
    width = 0.65,
    outlier.alpha = 0.35,
    colour = "grey25"
  ) +
  # The white diamond shows the mean NRMSE in each boxplot.
  ggplot2::stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    size = 2.3,
    fill = "white",
    colour = "black"
  ) +
  ggplot2::facet_grid(scenario ~ mechanism) +
  ggplot2::scale_fill_manual(
    values = c(
      "Hybrid kNN"   = "#E69F00",
      "Standard kNN" = "#CC3311",
      "Mean"         = "grey65"
    )
  ) +
  ggplot2::labs(
    x = NULL,
    y = "NRMSE",
    fill = NULL
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    legend.position = "none",
    strip.background = ggplot2::element_rect(fill = "grey92"),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    panel.grid.minor = ggplot2::element_blank()
  )

print(p_nrmse)

ggplot2::ggsave(
  filename = file.path(project_directory, "fig_nrmse_boxplots.pdf"),
  plot = p_nrmse,
  width = 9,
  height = 5.5
)

message("Created fig_nrmse_boxplots.pdf")
