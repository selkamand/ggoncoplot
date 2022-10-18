

# Oncoplot ----------------------------------------------------------------

#' GG oncoplot
#'
#' @param col_genes name of \strong{data} column containing gene names/symbols (string)
#' @param col_samples name of \strong{data} column containing sample identifiers (string)
#' @param col_mutation_type name of \strong{data} column describing mutation types (string)
#' @param col_tooltip name of \strong{data} column containing whatever information you want to display in (string)
#' @param topn how many of the top genes to visualise. Ignored if \code{genes_to_include} is supplied (number)
#' @param show_sample_ids show sample_ids_on_x_axis (flag)
#' @param .data data for oncoplot. A data.frame with 1 row per mutation in your cohort. Must contain columns describing gene_symbols and sample_identifiers, (data.frame)
#' @param genes_to_include specific genes to include in the oncoplot (character)
#' @param interactive should plot be interactive (boolean)
#' @param interactive_svg_width dimensions of interactive plot (number)
#' @param interactive_svg_height dimensions of interactive plot (number)
#' @param xlab_title x axis lable (string)
#' @param ylab_title y axis of interactive plot (number)
#' @param palette a named vector mapping all possible mutation types (vector names) to colours (vector values).
#' If not supplied ggoncoplot will check if all values are either valid SO or MAF variant classification terms
#' and use pre-made colour schemes for each of these ontologies from the \strong{mutationtypes} package.
#' If mutation type terms are not described using these ontologies, a 12 colour RColourBrewer palette will be used, but the user warned to make a custom mapping to force consistent colour schemes between plots (character)
#' @param fontsize_xlab size of x axis title (number)
#' @param fontsize_ylab size of y axis title (number)
#' @param fontsize_genes size of y axis text (gene names) (number)
#' @param fontsize_samples size of x axis text (sample names). Ignored unless show_sample_ids is set to true (number)
#' @return ggplot or girafe object if \code{interactive=TRUE}
#' @export
#'
#' @examples
#' # ===== GBM =====
#' gbm_csv <- system.file(
#'   package = "ggoncoplot",
#'   "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
#' )
#'
#' gbm_df <- read.csv(file = gbm_csv, header = TRUE)
#'
#' ggoncoplot(
#'   gbm_df,
#'   "Hugo_Symbol",
#'   "Tumor_Sample_Barcode",
#'   col_mutation_type = "Variant_Classification"
#' )
#'
ggoncoplot <- function(.data,
                       col_genes,
                       col_samples,
                       col_mutation_type = NULL,
                       genes_to_include = NULL,
                       col_tooltip = col_samples,
                       topn = 10,
                       palette = NULL,
                       show_sample_ids = FALSE,
                       interactive = TRUE,
                       interactive_svg_width = 12,
                       interactive_svg_height = 6,
                       xlab_title = "Sample",
                       ylab_title = "Gene",
                       fontsize_xlab = 26,
                       fontsize_ylab = 26,
                       fontsize_genes = 16,
                       fontsize_samples = 12) {
  assertthat::assert_that(is.data.frame(.data))
  assertthat::assert_that(assertthat::is.string(col_genes))
  assertthat::assert_that(assertthat::is.string(col_samples))
  assertthat::assert_that(is.null(col_mutation_type) | assertthat::is.string(col_mutation_type))
  assertthat::assert_that(is.null(genes_to_include) | is.character(genes_to_include))
  assertthat::assert_that(assertthat::is.string(col_tooltip))
  assertthat::assert_that(assertthat::is.number(topn))
  assertthat::assert_that(assertthat::is.number(fontsize_xlab))
  assertthat::assert_that(assertthat::is.number(fontsize_ylab))
  assertthat::assert_that(assertthat::is.number(fontsize_genes))
  assertthat::assert_that(assertthat::is.number(fontsize_samples))

  # Get dataframe with 1 row per sample-gene pair
  data_top_df <- ggoncoplot_prep_df(
    .data = .data,
    col_genes = col_genes, col_samples = col_samples,
    col_mutation_type = col_mutation_type,
    col_tooltip = col_tooltip,
    topn = topn,
    show_sample_ids = show_sample_ids,
    genes_to_include = genes_to_include
  )

  gg <- ggoncoplot_plot(
    .data = data_top_df,
    show_sample_ids = show_sample_ids, interactive = interactive, interactive_svg_width = interactive_svg_width,
    palette = palette,
    interactive_svg_height = interactive_svg_height,
    xlab_title = xlab_title,
    ylab_title = ylab_title,
    fontsize_xlab = fontsize_xlab,
    fontsize_ylab = fontsize_ylab,
    fontsize_genes = fontsize_genes,
    fontsize_samples = fontsize_samples
  )

  return(gg)
}


#' Prep data for oncoplot
#'
#' @param col_genes name of \strong{data} column containing gene names/symbols (string)
#' @param col_samples name of \strong{data} column containing sample identifiers (string)
#' @param col_mutation_type name of \strong{data} column describing mutation types (string)
#' @param col_tooltip name of \strong{data} column containing whatever information you want to display in (string)
#' @param topn how many of the top genes to visualise. If two genes are mutated in the same # of patients, 1 will be selected based on which appears first in the dataset.Ignored if \code{genes_to_include} is supplied (number)
#' @param show_sample_ids show sample_ids_on_x_axis (flag)
#' @param .data data for oncoplot. A data.frame with 1 row per mutation in your cohort. Must contain columns describing gene_symbols and sample_identifiers, (data.frame)
#' @param genes_to_include specific genes to include in the oncoplot (character)
#' @return dataframe with the following columns: 'Gene', 'Sample', 'MutationType', 'Tooltip'.
#' Sample is a factor with levels sorted in appropriate order for oncoplot vis.
#' Genes represents either topn genes or specific genes set by \code{genes_to_include}
#'
#' @examples
#' #' # ===== GBM =====
#' gbm_csv <- system.file(
#'   package = "ggoncoplot",
#'   "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
#' )
#'
#' gbm_df <- read.csv(file = gbm_csv, header = TRUE)
#'
#' ggoncoplot:::ggoncoplot_prep_df(
#'   gbm_df,
#'   col_genes = "Hugo_Symbol",
#'   col_samples = "Tumor_Sample_Barcode",
#'   col_mutation_type = "Variant_Classification"
#' )
#'
ggoncoplot_prep_df <- function(.data,
                               col_genes,
                               col_samples,
                               col_mutation_type = NULL,
                               col_tooltip = col_samples,
                               topn = 10,
                               show_sample_ids = FALSE,
                               genes_to_include = NULL) {
  assertthat::assert_that(is.data.frame(.data))
  assertthat::assert_that(assertthat::is.string(col_genes))
  assertthat::assert_that(assertthat::is.string(col_samples))
  assertthat::assert_that(is.null(col_mutation_type) | assertthat::is.string(col_mutation_type))
  assertthat::assert_that(is.null(genes_to_include) | is.character(genes_to_include))
  assertthat::assert_that(assertthat::is.string(col_tooltip))
  assertthat::assert_that(assertthat::is.number(topn))


  # Check specified columns are in .data
  data_colnames <- names(.data)

  check_valid_dataframe_column(
    data = .data,
    colnames = c(
      col_samples,
      col_genes,
      col_tooltip
    )
  )

  # Check optional columns are in .data
  if (!is.null(col_mutation_type)) {
    check_valid_dataframe_column(data = .data, colnames = col_mutation_type)
  }

  # Ensure Sample Column is A factor
  .data[[col_samples]] <- as.factor(.data[[col_samples]])

  # Look exclusively at a custom set of genes
  if (!is.null(genes_to_include)) {
    genes_not_found <- genes_to_include[!genes_to_include %in% .data[[col_genes]]]

    if (length(genes_not_found) > 0) {
      cli::cli_alert_warning("Failed to find the following [{length(genes_not_found)}] genes in your dataset")

      cli::cli_alert("{genes_not_found}")
      cli::cli_alert_warning("Either no samples have mutations in the above genes, or you've got the wrong gene names")
    }

    if (length(genes_not_found) == length(genes_to_include)) {
      cli::cli_abort("Couldn't find any of the genes you supplied in your dataset. Either no samples have mutations in these genes, or you've got the wrong gene names")
    }

    genes_in_selected_set <- .data[[col_genes]] %in% genes_to_include

    .data <- .data |>
      dplyr::filter(.data[[col_genes]] %in% genes_to_include)
  }

  # Filter out genes_to_ignore
  ## I'll have to implement later :)


  # Identify top genes by frequency
  data_gene_counts <- .data |>
    dplyr::distinct(.data[[col_samples]], .data[[col_genes]], .keep_all = TRUE) |> # This line stops multiple mutations of the same gene in the same sample counting multiple times towards the mutation frequency.
    dplyr::count(.data[[col_genes]])

  # If use hasn't specifically specified genes to include, filter for the 'topn' most mutated genes
  if (!is.null(genes_to_include)) {
    data_top_genes_df <- data_gene_counts
  } else {
    data_top_genes_df <- data_gene_counts |>
      dplyr::slice_max(.data$n, n = topn, with_ties = FALSE) # Set with_ties = TRUE to allow topn genes to return extra genes if there are ties in # of samples mutated
  }

  # Rank Genes based on how many samples they're mutated in
  data_top_genes <- data_top_genes_df |>
    dplyr::pull(.data[[col_genes]])

  data_top_genes_rank <- rank(data_top_genes_df[["n"]], ties.method = "first")

  # Filter dataset to only include the topn genes
  data_top_df <- .data |>
    dplyr::filter(.data[[col_genes]] %in% data_top_genes)

  # Order Genes Variable
  data_top_df[[col_genes]] <- forcats::fct_relevel(data_top_df[[col_genes]], data_top_genes[order(data_top_genes_rank)])

  # Sort Samples by mutated gene
  data_top_df <- data_top_df |>
    dplyr::group_by(.data[[col_samples]]) |>
    dplyr::mutate(
      SampleRankScore = score_based_on_gene_rank(mutated_genes = .data[[col_genes]], genes_informing_score = data_top_genes, gene_rank = data_top_genes_rank) # add secondary ranking based on secondary
    ) |>
    dplyr::ungroup()
  data_top_df[[col_samples]] <- forcats::fct_rev(forcats::fct_reorder(data_top_df[[col_samples]], data_top_df$SampleRankScore))

  # Consolidate to 1 row per sample-gene combo (collapse multiple mutations per gene into 1 row)
  # If col_mutation_type is supplied, will set mutation type to 'Multiple' for genes mutated multiple times in one patient
  if (!is.null(col_mutation_type)) {
    data_top_df <- data_top_df |> # Need to figure out how to fix this.
      dplyr::group_by(.data[[col_samples]], .data[[col_genes]]) |>
      dplyr::summarise(
        MutationType = dplyr::case_when(
          is.null(col_mutation_type) ~ NA_character_,
          dplyr::n_distinct(.data[[col_mutation_type]]) > 1 ~ "Multiple",
          TRUE ~ unique(.data[[col_mutation_type]])
        ) |> unique() |> paste0(collapse = "; "),
        Tooltip = paste0(unique(.data[[col_tooltip]]), collapse = "; ") # Edit this line to change how tooltips are collapsed
      ) |>
      dplyr::ungroup()
  } else {
    data_top_df <- data_top_df |> # Need to figure out how to fix this.
      dplyr::group_by(.data[[col_samples]], .data[[col_genes]]) |>
      dplyr::summarise(
        MutationType = NA_character_,
        Tooltip = paste0(unique(.data[[col_tooltip]]), collapse = "; ") # Edit this line to change how tooltips are collapsed
      ) |>
      dplyr::ungroup()
  }


  # Select just the columns we need,
  data_top_df <- data_top_df |>
    dplyr::select(Sample = {{ col_samples }}, Gene = {{ col_genes }}, MutationType = .data[["MutationType"]], Tooltip = .data[["Tooltip"]])

  return(data_top_df)
}



#' Plot oncoplot
#'
#' This function takes the output from \strong{ggoncoplot_prep_df} and plots it.
#' Should not be exposed since it makes some assumptions about structure of input data.
#'
#' @inheritParams ggoncoplot
#' @param .data transformed data from [ggoncoplot_prep_df()] (data.frame)
#' @inherit ggoncoplot return
#' @inherit ggoncoplot examples
ggoncoplot_plot <- function(.data,
                            show_sample_ids = FALSE,
                            interactive = TRUE,
                            palette = NULL,
                            interactive_svg_width = 12,
                            interactive_svg_height = 6,
                            xlab_title = "Sample",
                            ylab_title = "Gene",
                            fontsize_xlab = 16,
                            fontsize_ylab = 16,
                            fontsize_genes = 14,
                            fontsize_samples = 10
                            ) {
  check_valid_dataframe_column(.data, c("Gene", "Sample", "MutationType", "Tooltip"))

  # Consistent Colour Scheme
  unique_impacts <- unique(.data[["MutationType"]])
  unique_impacts_minus_multiple <- unique_impacts[unique_impacts != "Multiple"]


  if (all(is.na(unique_impacts))) {
    palette <- NA
  } else if (is.null(palette)) {
    mutation_dictionary <- mutationtypes::mutation_types_identify(unique_impacts_minus_multiple)

    if (mutation_dictionary == "MAF") {
      cli::cli_alert_success("Mutation Types are described using valid MAF terms ... using MAF palete")
      palette <- c(mutationtypes::mutation_types_maf_palette(), Multiple = "black")
      palette <- palette[names(palette) %in% unique_impacts]
    } else if (mutation_dictionary == "SO") {
      cli::cli_alert_success("Mutation Types are described using valid SO terminology ... using SO palete")
      palette <- c(mutationtypes::mutation_types_so_palette(), Multiple = "black")
      palette <- palette[names(palette) %in% unique_impacts]
    } else {
      cli::cli_alert_warning("Mutation Types are not described with any known ontology.
                             Using an RColorBrewer palette by default.
                             When running this plot with other datasets, it is possible the colour scheme may differ.
                             We STRONGLY reccomend supplying a custom MutationType -> colour mapping using the {.arg palette} argument")

      # .data[['MutationType']] <- forcats::fct_infreq(f = .data[['MutationType']])
      rlang::check_installed("RColorBrewer", reason = "To create default palette for `ggoncoplot()`")
      palette <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
    }
  } else { # What if custom palette is supplied?
    if (!all(unique_impacts %in% names(palette))) {
      terms_without_mapping <- unique_impacts[!unique_impacts %in% names(palette)]
      cli::cli_abort("Please add colour mappings for the following terms: {terms_without_mapping}")
      palette <- palette[names(palette) %in% unique_impacts]
    }
  }


  # Create ggplot
  gg <- ggplot2::ggplot(
    data = .data,
    mapping = ggplot2::aes_string(
      y = "Gene",
      x = "Sample",
      fill = "MutationType"
    )
  )

  # Add interactive/non-interactive geom layer
  #browser()
  # Todo: rework this to avoid tidyr dependency. Maybe .data shouldl have an entry for each level.
  #browser()
  cli::cli_alert_info('')
  #browser()
  gg <- gg +
    ggiraph::geom_tile_interactive(
      data=expand.grid(Sample = droplevels(unique(.data[['Sample']])), Gene = unique(.data[["Gene"]])),
      ggplot2::aes_string(
        tooltip = "Sample", # Can't just use tooltip since these don't have a value in .data. Maybe I should fix the source problem
        data_id = "Sample",
      ),
      height = 0.9,
      width = 0.95,
      fill = "grey90"
    ) +
    ggiraph::geom_tile_interactive(
      data = .data,
      ggplot2::aes_string(
        tooltip = "Tooltip",
        data_id = "Sample",
        height = 0.9,
        width = 0.95
      )
    )



  # Label axis
  gg <- gg + ggplot2::xlab(xlab_title) + ggplot2::ylab(ylab_title)

  # Add fill colour
  gg <- gg + ggplot2::scale_fill_manual(values = palette)


  # Apply default theme
  gg <- gg + theme_oncoplot_default()

  # Add line between genes
  gg <- gg + ggplot2::geom_hline(yintercept = seq(0, length(unique(.data[['Gene']]))) + .5, color="gray30")

  # Change text size for x and y axis labels
  gg <- gg + ggplot2::theme(
    axis.title.x = ggplot2::element_text(size = fontsize_xlab),
    axis.title.y = ggplot2::element_text(size = fontsize_ylab),
    axis.text.x  = ggplot2::element_text(size = fontsize_samples, angle = 45, hjust = 1),
    axis.text.y  = ggplot2::element_text(size = fontsize_genes),
    axis.title = ggplot2::element_text(face = "bold")
  )

  # Panel changes
  gg <- gg + ggplot2::theme(
    panel.grid.major = ggplot2::element_blank()
  )


  # Show/hide sample ids on x axis
  if (!show_sample_ids) {
    gg <- gg + ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    )
  }

  # Turn gg into an interactive ggiraph object if interactive = TRUE
  if (interactive) {
    gg <- ggiraph::girafe(
      width_svg = interactive_svg_width, height_svg = interactive_svg_height,
      ggobj = gg,
      options = list(
        ggiraph::opts_tooltip(
          opacity = .8,
          css = "background-color:gray;color:white;padding:2px;border-radius:2px;"
        ),
        ggiraph::opts_hover_inv(css = "opacity:0.2;"),
        ggiraph::opts_hover(css = "stroke-width:5;"),
        ggiraph::opts_selection("stroke-width:5;opacity:1", type = "multiple", only_shiny = FALSE)

      )
    )
  }
  return(gg)
}


# Utils -------------------------------------------------------------------

#' Generate score based on genes
#'
#' Score used to sort samples based on which genes are mutated. Make sure to run on one sample at once (use grouping)
#'
#' @param mutated_genes vector of genes that are mutated for a single sample (character)
#' @param genes_informing_score which genes determine the sort order? (character)
#' @param gene_rank what is the order of importance of genes used to determine sort order. Higher number = higher in sort order (character)
#' @param debug_mode debug mode (flag)
#'
#' @return a score (higher = should be higher in the sorting order) (number)
#'
#' @examples
#' \dontrun{
#' # First set of genes has a high rank since both BRCA2 and EGFR are mutated
#' score_based_on_gene_rank(c("TERT", "EGFR", "PTEN", "BRCA2"), c("EGFR", "BRCA2"), gene_rank = 1:2)
#'
#' # If EGFR is mutated without BRCA2, we get a lower score
#' score_based_on_gene_rank(c("TERT", "EGFR", "PTEN", "IDH1"), c("EGFR", "BRCA2"), gene_rank = 1:2)
#'
#' # If BRCA2 is mutated without EGFR,
#' # we get a score lower than BRCA2+EGFR but higher than EGFR alone due to higher gene_rank of BRCA2
#' score_based_on_gene_rank(c("TERT", "IDH1", "PTEN", "BRCA2"), c("EGFR", "BRCA2"), gene_rank = 1:2)
#' }
score_based_on_gene_rank <- function(mutated_genes, genes_informing_score, gene_rank, debug_mode = FALSE) {
  assertthat::assert_that(is.character(mutated_genes) | is.factor(mutated_genes))
  assertthat::assert_that(is.character(genes_informing_score))
  assertthat::assert_that(is.numeric(gene_rank))
  assertthat::assert_that(length(genes_informing_score) == length(gene_rank))


  gene_order <- rank(gene_rank, ties.method = "first")


  res <- genes_informing_score %in% mutated_genes
  names(res) <- genes_informing_score
  base10_values <- ifelse(res, yes = 2^(gene_order - 1), no = 0)

  total_score <- sum(base10_values)

  if (debug_mode) {
    cli::cli_alert("base_values = [{names(res)}].")
    cli::cli_alert("base_values = [{base10_values}].")
    cli::cli_alert("total score = [{total_score}].")
  }
  return(total_score)
}


#' Oncoplot Theme: default
#'
#' @param ... passed to [ggplot2::theme()] theme
#' @importFrom ggplot2 %+replace%
theme_oncoplot_default <- function(...) {
  ggplot2::theme_bw(...) %+replace%
    ggplot2::theme(
      panel.border = ggplot2::element_rect(size = 1, fill = NA),
      #panel.grid.minor = ggplot2::element_line(colour = "red"),
      panel.grid.major = ggplot2::element_blank(),
      # panel.grid.minor.y = ggplot2::element_line(colour = "red"),
      axis.title = ggplot2::element_text(face = "bold")
    )
}


#' data.frame has colnames
#'
#' Assert that data.frame contains a set of user defined column names.
#'
#' data.frame may have any additional colnames.
#' It just has to have AT LEAST the columns specified by \code{colnames}
#'
#' @param data dataframe that you want to assert contain specific columns (data.frame)
#' @param colnames Name (character)
#' @param error_call error call environment (do not change)
#'
#' @return invisibly returns TRUE. If data is missing columns, will throw error
#'
#' @examples
#' # Check mtcars has columns 'mpg' and 'cyl'
#' ggoncoplot:::check_valid_dataframe_column(mtcars, c("mpg", "cyl"))
#'
#' @details Informs user about the missing columns one at a time. This may change in future
#'
check_valid_dataframe_column <- function(data, colnames, error_call = rlang::caller_env()) {
  assertthat::assert_that(is.character(colnames))
  assertthat::assert_that(is.data.frame(data))


  data_colnames <- colnames(data)

  for (colname in colnames) {
    if (!colname %in% data_colnames) {
      cli::cli_abort(
        c(
          "!" = "Could find column: [{colname}] in input data",
          "!" = "Please select one of [{data_colnames}]"
        )
      )
    }
  }
  invisible(TRUE)
}
