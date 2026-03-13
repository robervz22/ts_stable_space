remove(list = ls())

# ============================================================
# Empirical illustrations for:
# "On the Stable Space of a Multivariate Time Series"
#
# This script generates, displays, and exports:
# 1) Observed vs fitted values by method
# 2) Projection accuracy metrics by method
#
# Output formats:
# - EPS (for LaTeX / journal submission)
# - PDF (for vector output and sharing)
# - PNG (for quick inspection)
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(grid)

# ============================================================
# User options
# ============================================================

#################
# ggplot2 setup #
#################
mysize <- 14 # base font size
mytheme <- theme_bw(base_size = mysize) + 
  theme(
    axis.title = element_text(size = rel(1.2), face = "bold"),
    axis.text = element_text(size = rel(1), face = "bold"),
    legend.title = element_text(size = rel(1.1), face = "bold"),
    legend.text = element_text(size = rel(1)),
    plot.title = element_text(size = rel(1.4), hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(size = rel(1.2), hjust = 0.5, face = "bold"),
    plot.caption = element_text(size = rel(0.8), hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )

ggadd <- function(plot, object) {
  update_ggplot(object, plot)
}


path <- getwd()

show_x11 <- FALSE
show_covid_line <- TRUE
covid_date <- as.Date("2020-03-01")

# TRUE  -> cumulative RMSE: sqrt( cumsum((obs-fit)^2) / t )
# FALSE -> pointwise scaled error: sqrt( (obs-fit)^2 / H )
use_cumulative_rmse <- TRUE

eps_width  <- 14
eps_height <- 7
eps_family <- "Times"

png_width  <- 14
png_height <- 7
png_dpi    <- 300

# Data
inflation <- read.csv(paste(path, "/databases/fit_inflation.csv", sep = ""), row.names = 1)
bie <- read.csv(paste(path, "/databases/fit_BIE.csv", sep = ""), row.names = 1)

# ============================================================
# Main reusable function
# ============================================================

build_projection_plots <- function(
    data,
    series_prefix,
    path,
    output_tag,
    show_x11 = TRUE,
    show_covid_line = TRUE,
    covid_date = as.Date("2020-03-01"),
    use_cumulative_rmse = TRUE,
    eps_width = 14,
    eps_height = 7,
    eps_family = "Times",
    png_width = 14,
    png_height = 7,
    png_dpi = 300
) {
  
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  
  # ----------------------------------------------------------
  # 1) Prepare data
  # ----------------------------------------------------------
  
  df <- data %>%
    rownames_to_column(var = "Date") %>%
    mutate(Date = as.Date(Date))
  
  H <- nrow(df)
  
  obs_cols <- grep("_obs$", names(df), value = TRUE)
  
  obs_cols <- obs_cols[sapply(obs_cols, function(x) {
    any(startsWith(x, paste0(series_prefix, "_")))
  })]
  
  if (length(obs_cols) == 0) {
    stop("No observed columns matching series_prefix were found.")
  }
  
  variables <- sub("_obs$", "", obs_cols)
  
  fitted_cols_all <- unlist(lapply(variables, function(v) {
    grep(paste0("^", v, "_"), names(df), value = TRUE)
  }), use.names = FALSE)
  
  fitted_cols_all <- setdiff(fitted_cols_all, obs_cols)
  detected_methods <- unique(sub("^.*?_", "", fitted_cols_all))
  
  preferred_methods <- c(
    "stationary_PLS",
    "stationary_PCA",
    "stationary_SPLS",
    "stationary_SPCA",
    "classic_PLS",
    "classic_PCA",
    "Johansen"
  )
  
  methods <- preferred_methods[preferred_methods %in% detected_methods]
  
  if (length(methods) == 0) {
    stop("No fitted methods were detected.")
  }
  
  method_labels_map <- c(
    stationary_PLS  = "Stationary PLS",
    stationary_PCA  = "Stationary PCA",
    stationary_SPLS = "Stationary SPLS",
    stationary_SPCA = "Stationary SPCA",
    classic_PLS     = "Classic PLS",
    classic_PCA     = "Classic PCA",
    Johansen        = "Johansen"
  )
  
  method_labels <- unname(method_labels_map[methods])
  
  # ----------------------------------------------------------
  # 2) Long data
  # ----------------------------------------------------------
  
  long_df <- bind_rows(
    lapply(variables, function(v) {
      obs_col <- paste0(v, "_obs")
      fit_cols <- paste0(v, "_", methods)
      fit_cols <- fit_cols[fit_cols %in% names(df)]
      
      df %>%
        select(Date, all_of(obs_col), all_of(fit_cols)) %>%
        pivot_longer(
          cols = all_of(fit_cols),
          names_to = "series",
          values_to = "fitted"
        ) %>%
        mutate(
          variable = v,
          method   = sub(paste0("^", v, "_"), "", series),
          observed = .data[[obs_col]],
          error    = observed - fitted
        ) %>%
        select(Date, variable, method, observed, fitted, error)
    })
  ) %>%
    mutate(
      variable = factor(variable, levels = variables),
      method   = factor(method, levels = methods, labels = method_labels)
    ) %>%
    group_by(variable, method) %>%
    arrange(Date, .by_group = TRUE) %>%
    mutate(
      t = row_number(),
      metric = if (use_cumulative_rmse) {
        sqrt(cumsum(error^2) / t)
      } else {
        sqrt((error^2) / H)
      }
    ) %>%
    ungroup()
  
  # ----------------------------------------------------------
  # 3) Plot A: observed vs fitted
  # ----------------------------------------------------------
  
  p_fit <- ggplot(long_df, aes(x = Date))
  p_fit <- ggadd(p_fit, geom_line(aes(y = observed), linewidth = 0.60))
  p_fit <- ggadd(p_fit, geom_line(aes(y = fitted), linewidth = 0.55, linetype = 2))
  p_fit <- ggadd(p_fit, facet_grid(variable ~ method, scales = "free_y"))
  p_fit <- ggadd(p_fit, labs(
    x = "Date",
    y = NULL
  ))
  p_fit <- ggadd(p_fit, mytheme)
  p_fit <- ggadd(p_fit, theme(
    # strip.background = element_rect(fill = "white", colour = "black"),
    # panel.grid.minor = element_blank(),
    # panel.spacing = unit(0.55, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # plot.margin = margin(8, 8, 8, 8)
  ))
  
  if (show_covid_line) {
    p_fit <- ggadd(p_fit, geom_vline(
      xintercept = covid_date,
      linetype = 3,
      linewidth = 0.4
    ))
  }
  
  # ----------------------------------------------------------
  # 4) Plot B: metric
  # ----------------------------------------------------------
  
  y_lab_metric <- if (use_cumulative_rmse) {
    "Cumulative RMSE"
  } else {
    expression(sqrt((obs - hat(obs))^2 / H))
  }
  
  p_metric <- ggplot(long_df, aes(x = Date, y = metric))
  p_metric <- ggadd(p_metric, geom_line(linewidth = 0.55))
  p_metric <- ggadd(p_metric, facet_grid(variable ~ method, scales = "free_y"))
  p_metric <- ggadd(p_metric, labs(
    x = "Date",
    y = y_lab_metric
  ))
  p_metric <- ggadd(p_metric, mytheme)
  p_metric <- ggadd(p_metric, theme(
    # strip.background = element_rect(fill = "white", colour = "black"),
    # panel.grid.minor = element_blank(),
    # panel.spacing = unit(0.55, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    # plot.margin = margin(8, 8, 8, 8)
  ))
  
  if (show_covid_line) {
    p_metric <- ggadd(p_metric, geom_vline(
      xintercept = covid_date,
      linetype = 3,
      linewidth = 0.4
    ))
  }
  
  # ----------------------------------------------------------
  # 5) Show on screen
  # ----------------------------------------------------------
  
  if (show_x11) {
    x11(width = eps_width, height = eps_height)
    print(p_fit)
    
    x11(width = eps_width, height = eps_height)
    print(p_metric)
  }
  
  # ----------------------------------------------------------
  # 6) File names
  # ----------------------------------------------------------
  
#   fit_eps_file <- file.path(
#     path,
#     paste0("projection_observed_vs_fitted_", output_tag, ".eps")
#   )
  
#   fit_png_file <- file.path(
#     path,
#     paste0("projection_observed_vs_fitted_", output_tag, ".png")
#   )
  
#   fit_pdf_file <- file.path(
#     path,
#     paste0("projection_observed_vs_fitted_", output_tag, ".pdf")
#   )
  
  metric_eps_file <- file.path(
    path,
    if (use_cumulative_rmse) {
      paste0("images/Fig", output_tag, ".eps")
    } else {
      paste0("Fig", output_tag, ".eps")
    }
  )
  
  metric_png_file <- file.path(
    path,
    if (use_cumulative_rmse) {
      paste0("images/Fig", output_tag, ".png")
    } else {
      paste0("images/Fig", output_tag, ".png")
    }
  )
  
  metric_pdf_file <- file.path(
    path,
    if (use_cumulative_rmse) {
      paste0("images/Fig", output_tag, ".pdf")
    } else {
      paste0("images/Fig", output_tag, ".pdf")
    }
  )
  
  # ----------------------------------------------------------
  # 7) Export EPS
  # ----------------------------------------------------------
  
#   ggsave(
#     filename = fit_eps_file,
#     plot = p_fit,
#     device = cairo_ps,
#     fallback_resolution = 600,
#     width = eps_width,
#     height = eps_height,
#     units = "in",
#     family = eps_family
#   )
  
  ggsave(
    filename = metric_eps_file,
    plot = p_metric,
    device = cairo_ps,
    fallback_resolution = 600,
    width = eps_width,
    height = eps_height,
    units = "in"
    # family = eps_family
  )
  
  # ----------------------------------------------------------
  # 8) Export PDF
  # ----------------------------------------------------------
  
#   ggsave(
#     filename = fit_pdf_file,
#     plot = p_fit,
#     device = cairo_pdf,
#     width = eps_width,
#     height = eps_height,
#     units = "in",
#     family = eps_family
#   )
  
  ggsave(
    filename = metric_pdf_file,
    plot = p_metric,
    device = cairo_pdf,
    width = eps_width,
    height = eps_height,
    units = "in"
    # family = eps_family
  )
  
  # ----------------------------------------------------------
  # 9) Export PNG
  # ----------------------------------------------------------
  
#   ggsave(
#     filename = fit_png_file,
#     plot = p_fit,
#     width = png_width,
#     height = png_height,
#     units = "in",
#     dpi = png_dpi
#   )
  
  ggsave(
    filename = metric_png_file,
    plot = p_metric,
    width = png_width,
    height = png_height,
    units = "in",
    dpi = png_dpi
  )
  
  # ----------------------------------------------------------
  # 10) Messages
  # ----------------------------------------------------------
  
#   message("Observed vs fitted EPS exported to: ", fit_eps_file)
#   message("Observed vs fitted PDF exported to: ", fit_pdf_file)
#   message("Observed vs fitted PNG exported to: ", fit_png_file)
  message("Metric EPS exported to: ", metric_eps_file)
  message("Metric PDF exported to: ", metric_pdf_file)
  message("Metric PNG exported to: ", metric_png_file)
  
  invisible(list(
    observed_vs_fitted = p_fit,
    metric_plot = p_metric,
    data_long = long_df,
    methods_detected = methods,
    variables_detected = variables,
    exported_files = c(
    #   fit_eps_file,
    #   fit_pdf_file,
    #   fit_png_file,
      metric_eps_file,
      metric_pdf_file,
      metric_png_file
    )
  ))
}

# ============================================================
# Run for inflation
# ============================================================

plots_inflation <- build_projection_plots(
  data = inflation,
  series_prefix = c("P", "W"),
  path = path,
  output_tag = "6",
  show_x11 = show_x11,
  show_covid_line = show_covid_line,
  covid_date = covid_date,
  use_cumulative_rmse = use_cumulative_rmse,
  eps_width = eps_width,
  eps_height = eps_height,
  eps_family = eps_family,
  png_width = png_width,
  png_height = png_height,
  png_dpi = png_dpi
)

# ============================================================
# Run for BIE
# ============================================================

plots_bie <- build_projection_plots(
  data = bie,
  series_prefix = c("s23", "s54"),
  path = path,
  output_tag = "7",
  show_x11 = show_x11,
  show_covid_line = show_covid_line,
  covid_date = covid_date,
  use_cumulative_rmse = use_cumulative_rmse,
  eps_width = eps_width,
  eps_height = eps_height,
  eps_family = eps_family,
  png_width = png_width,
  png_height = png_height,
  png_dpi = png_dpi
)
