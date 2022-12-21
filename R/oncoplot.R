

# Oncoplot ----------------------------------------------------------------

#' GG oncoplot
#' @importFrom patchwork plot_layout
#'
#' @param col_genes name of **data** column containing gene names/symbols (string)
#' @param col_samples name of **data** column containing sample identifiers (string)
#' @param col_mutation_type name of **data** column describing mutation types (string)
#' @param col_tooltip name of **data** column containing whatever information you want to display in (string)
#' @param topn how many of the top genes to visualise. Ignored if `genes_to_include` is supplied (number)
#' @param show_sample_ids show sample_ids_on_x_axis (flag)
#' @param .data data for oncoplot. A data.frame with 1 row per mutation in your cohort. Must contain columns describing gene_symbols and sample_identifiers (data.frame)
#' @param genes_to_include specific genes to include in the oncoplot (character)
#' @param genes_to_ignore names of the genes that should be ignored (character)
#' @param return_extra_genes_if_tied instead of strictly returning `topn` genes,
#' in the case of ties (where multiple genes are mutated in the exact same number of samples, complicating selection of top n genes), return all tied genes (potentially more than topn).
#' If FALSE, will return strictly `topn` genes, breaking ties based on order of appearance in dataset (flag)
#' @param interactive should plot be interactive (boolean)
#' @param interactive_svg_width dimensions of interactive plot (number)
#' @param interactive_svg_height dimensions of interactive plot (number)
#' @param xlab_title x axis lable. Set `xlab_title = NULL` to remove title (string)
#' @param ylab_title y axis of interactive plot. Set `ylab_title = NULL` to remove title (string)
#' @param palette a named vector mapping all possible mutation types (vector names) to colours (vector values).
#' If not supplied ggoncoplot will check if all values are either valid SO or MAF variant classification terms
#' and use pre-made colour schemes for each of these ontologies from the **mutationtypes** package.
#' If mutation type terms are not described using these ontologies, a 12 colour RColourBrewer palette will be used, but the user warned to make a custom mapping to force consistent colour schemes between plots (character)
#' @param fontsize_xlab size of x axis title (number)
#' @param fontsize_ylab size of y axis title (number)
#' @param fontsize_genes size of y axis text (gene names) (number)
#' @param fontsize_samples size of x axis text (sample names). Ignored unless show_sample_ids is set to true (number)
#' @param fontsize_tmb_title fontsize of y axis title for TMB marginal plot (number)
#' @param fontsize_tmb_axis fontsize of y axis text for TMB marginal plot (number)
#' @param verbose verbose mode (flag)
#' @param tile_height  proportion of available vertical space each tile will take up (0-1) (number)
#' @param tile_width proportion of available horizontal space  each tile take up (0-1) (number)
#' @param colour_backround colour used for background non-mutated tiles (string)
#' @param fontsize_count fontsize of gene mutation count x axis (number)
#'
#' @param draw_gene_barplot add a barplot describing number of samples with each gene mutated (right side) (flag)
#' @param draw_tmb_barplot add a barplot describing total number of mutations in each sample (above main plot). If a single gene is mutated multiple times, all mutations are counted towards total (flag)
#' @param show_all_samples show all samples in oncoplot,
#' @param log10_transform_tmb log10 transform total number of mutations for TMB marginal plot (flag)
#' even if they don't have mutations in the selected genes.
#' Samples only described in metadata but with no mutations
#' at all are still filtered out by default, but you can show these too by setting `metadata_require_mutations = FALSE`  (flag)
#' @param metadata_require_mutations filter out samples from metadata lacking any mutations in .data (flag)
#' @param metadata dataframe describing sample level metadata.
#' One column must contain unique sample identifiers. Other columns can describe numeric / categorical metadata (data.frame)
#' @param col_samples_metadata which column in metadata data.frame describes sample identifiers (string)
#' @param show_ylab_title show y axis title of oncoplot (flag)
#' @param show_xlab_title show x axis title of oncoplot (flag)
#' @param show_ylab_title_tmb show y axis title of TMB margin plot (flag)
#' @param colour_mutation_type_unspecified colour of mutations in oncoplot and margin plots if `col_mutation_type` is not supplied (string)
#' @return ggplot or girafe object if `interactive=TRUE`
#' @export
#'
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
                       genes_to_ignore = NULL,
                       col_tooltip = col_samples,
                       topn = 10,
                       return_extra_genes_if_tied = FALSE,
                       palette = NULL,
                       metadata = NULL,
                       col_samples_metadata = col_samples,
                       metadata_require_mutations = TRUE,
                       show_sample_ids = FALSE,
                       interactive = TRUE,
                       interactive_svg_width = 12,
                       interactive_svg_height = 6,
                       xlab_title = "Sample",
                       ylab_title = "Gene",
                       fontsize_xlab = 26,
                       fontsize_ylab = 26,
                       fontsize_genes = 16,
                       fontsize_samples = 12,
                       fontsize_count = 14,
                       fontsize_tmb_title = 14,
                       fontsize_tmb_axis = 11,
                       tile_height = 1,
                       tile_width = 1,
                       colour_backround = "grey90",
                       colour_mutation_type_unspecified = "grey10",
                       draw_gene_barplot = FALSE,
                       draw_tmb_barplot = FALSE,
                       show_ylab_title = FALSE,
                       show_xlab_title = FALSE,
                       show_ylab_title_tmb = FALSE,
                       log10_transform_tmb = TRUE,
                       show_all_samples = FALSE,
                       verbose = TRUE
                       ) {


  # Assertions --------------------------------------------------------------
  assertthat::assert_that(is.data.frame(.data))
  assertthat::assert_that(nrow(.data) > 0)
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
  assertthat::assert_that(assertthat::is.number(fontsize_tmb_title))
  assertthat::assert_that(assertthat::is.number(fontsize_tmb_axis))
  assertthat::assert_that(assertthat::is.flag(verbose))
  assertthat::assert_that(assertthat::is.number(tile_height))
  assertthat::assert_that(assertthat::is.number(tile_width))
  assertthat::assert_that(assertthat::is.string(colour_backround))
  assertthat::assert_that(assertthat::is.flag(draw_gene_barplot))
  assertthat::assert_that(assertthat::is.flag(draw_tmb_barplot))
  assertthat::assert_that(is.null(metadata) | is.data.frame(metadata))
  assertthat::assert_that(!anyDuplicated(metadata[[col_samples_metadata]]))
  assertthat::assert_that(assertthat::is.flag(log10_transform_tmb))
  assertthat::assert_that(assertthat::is.string(colour_mutation_type_unspecified))


  # Configuration -----------------------------------------------------------
  # Properties we might want to tinker with, but not expose to user

  # Plot margins for tile plot
  # margins on right and top will be forced to zero if
  # marginal plots (TMB / gene barplots) are added
  margin_main_t = 0.2
  margin_main_r = 0.3
  margin_main_b = 0.2
  margin_main_l = 0.3
  margin_units = "pt"


  # Metadata preprocessing --------------------------------------------------

  # Remove any samples with metadata but ZERO mutations (can turn this off)
  if(metadata_require_mutations & !is.null(metadata)){
    lgl_samples_have_muts <- metadata[[col_samples_metadata]] %in% unique(.data[[col_samples]])
    samples_without_muts <- unique(metadata[[col_samples_metadata]][!lgl_samples_have_muts])

    if(verbose){
      cli::cli_alert_info("{length(samples_without_muts)} samples with metadata have no mutations. Fitering these out")
      cli::cli_alert_info("To keep these samples, set {.arg metadata_require_mutations = FALSE}. To view them in the oncoplot ensure you additionally set {.arg show_all_samples = TRUE}")
      names(samples_without_muts) <- rep(">", times = length(samples_without_muts))
      cli::cli_bullets(samples_without_muts)
    }
    metadata <- metadata[lgl_samples_have_muts,]
  }


  # Gene Order  --------------------------------------------------------------
  # Get Genes in Order for Oncoplot
  genes_for_oncoplot <- get_genes_for_oncoplot(
    .data = .data,
    col_samples = col_samples,
    col_genes = col_genes,
    topn = topn,
    genes_to_ignore = genes_to_ignore,
    return_extra_genes_if_tied = return_extra_genes_if_tied,
    genes_to_include = genes_to_include,
    verbose = verbose
  )



  # Preprocess dataframe ----------------------------------------------------
  # Get dataframe with 1 row per sample-gene pair
  data_top_df <- ggoncoplot_prep_df( # Add a samples_for_oncoplot
    .data = .data,
    col_genes = col_genes, col_samples = col_samples,
    col_mutation_type = col_mutation_type,
    col_tooltip = col_tooltip,
    genes_for_oncoplot = genes_for_oncoplot,
    verbose=verbose
  )


  # Sample order ----------------------------------------------
  # Get Sample Order,
  samples_with_mutations_in_selected_genes <- levels(droplevels(data_top_df[["Sample"]]))
  samples_with_any_mutations <- unique(.data[[col_samples]])

  # note we've already filtered out samples lacking any mutations above
  # (unless metadata_require_mutations == TRUE)
  samples_with_clinical_metadata <- metadata[[col_samples_metadata]]


  # The order of samples on x axis is determined by order in all_sample_ids
  # By default we keep order from `data_top_df` (mutation based ranking)
  # then tack on samples with mutations in unselected genes
  # then add samples with clinical metadata but no mutations
  all_sample_ids <- unique(c(
    samples_with_mutations_in_selected_genes,
    samples_with_any_mutations,
    samples_with_clinical_metadata
  ))

  if(!show_all_samples)
    samples_to_show <- samples_with_mutations_in_selected_genes
  else
    samples_to_show <- all_sample_ids

  # Add code for changing order of samples here
  # Example all_sample_ids = reorder_by_clinical_property(all_sample_ids, clinical_property)

  # Here we take each dataframe, ensure content only describes samples_to_show,
  # and any missing samples are added as factor levels.
  # This lets us just use scale_x_discrete(drop=FALSE) when plotting to show all samples we care about
  .data <- unify_samples(.data = .data, col_samples = col_samples, samples_to_show = samples_to_show)
  data_top_df <- unify_samples(.data = data_top_df, col_samples = "Sample", samples_to_show = samples_to_show)
  metadata <- unify_samples(.data = metadata, col_samples = col_samples_metadata, samples_to_show = samples_to_show)


  # Palette -----------------------------------------------------------------
  palette <- topn_to_palette(.data = data_top_df, palette = palette, verbose = verbose)
  #browser()


  # Draw main plot --------------------------------------------------------
  gg_main <- ggoncoplot_plot(
    .data = data_top_df,
    show_sample_ids = show_sample_ids,
    palette = palette,
    xlab_title = xlab_title,
    ylab_title = ylab_title,
    fontsize_xlab = fontsize_xlab,
    fontsize_ylab = fontsize_ylab,
    fontsize_genes = fontsize_genes,
    fontsize_samples = fontsize_samples,
    tile_height = tile_height,
    tile_width = tile_width,
    colour_backround = colour_backround,
    margin_t = margin_main_t,
    margin_r = margin_main_r,
    margin_b = margin_main_b,
    margin_l = margin_main_l,
    show_ylab_title = show_ylab_title,
    show_xlab_title = show_ylab_title,
    margin_unit = margin_units,
    colour_mutation_type_unspecified = colour_mutation_type_unspecified
  )



  # Draw marginal plots -----------------------------------------------------
  gg_gene_barplot = NULL
  gg_tmb_barplot = NULL

  ## Gene Barplot -----------------------------------------------------------
  if(draw_gene_barplot){
    gg_gene_barplot <- ggoncoplot_gene_barplot(
      .data = data_top_df,
      fontsize_count = fontsize_count,
      palette = palette,
      colour_mutation_type_unspecified = colour_mutation_type_unspecified
    )

  }

  ## TMB plot  -----------------------------------------------------------
  if(draw_tmb_barplot){
    gg_tmb_barplot <- ggoncoplot_tmb_barplot(
      .data = .data,
      col_samples = col_samples,
      col_mutation_type = col_mutation_type,
      log10_transform = log10_transform_tmb,
      fontsize_ylab = fontsize_tmb_title,
      fontsize_axis_text = fontsize_tmb_axis,
      show_ylab = show_ylab_title_tmb,
      palette = palette,
      colour_mutation_type_unspecified = colour_mutation_type_unspecified
    )

  }

  ## Combine marginal plots -----------------------------------------------------------
  gg_final <- combine_plots(gg_main, gg_tmb = gg_tmb_barplot, gg_gene = gg_gene_barplot)


  ## Draw metadata tiles ---------------------------------------------------------


  # Make Interactive -------------------------------------------------------
  # Turn gg into an interactive ggiraph object if interactive = TRUE
  if (interactive) {
    gg_final <- ggiraph::girafe(
      width_svg = interactive_svg_width, height_svg = interactive_svg_height,
      ggobj = gg_final,
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

  return(gg_final)
}



# Data Transformation -----------------------------------------------------


#' Prep data for oncoplot
#'
#' @inheritParams ggoncoplot
#' @param col_genes name of **data** column containing gene names/symbols (string)
#' @param col_samples name of **data** column containing sample identifiers (string)
#' @param col_mutation_type name of **data** column describing mutation types (string)
#' @param col_tooltip name of **data** column containing whatever information you want to display in (string)
#' @param .data data for oncoplot. A data.frame with 1 row per mutation in your cohort. Must contain columns describing gene_symbols and sample_identifiers (data.frame)
#' @param genes_for_oncoplot a list of genes to include in the oncoplot (character).
#' @return dataframe with the following columns: 'Gene', 'Sample', 'MutationType', 'Tooltip'.
#' Sample is a factor with levels sorted in appropriate order for oncoplot vis.
#' Genes represents either topn genes or specific genes set by `genes_to_include`
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
#' # Get genes in appropriate order for oncoplot
#' genes_for_oncoplot <- ggoncoplot:::get_genes_for_oncoplot(
#'   .data = gbm_df,
#'   col_samples = "Tumor_Sample_Barcode",
#'   col_genes = "Hugo_Symbol",
#'   topn = 20,
#'   verbose = FALSE
#' )
#'
#' # Create dataframe basis of oncoplot (1 row per sample-gene combo)
#' ggoncoplot:::ggoncoplot_prep_df(
#'   gbm_df,
#'   col_genes = "Hugo_Symbol",
#'   col_samples = "Tumor_Sample_Barcode",
#'   col_mutation_type = "Variant_Classification",
#'   genes_for_oncoplot = genes_for_oncoplot
#' )
#'
ggoncoplot_prep_df <- function(.data,
                               col_genes,
                               col_samples,
                               genes_for_oncoplot,
                               col_mutation_type = NULL,
                               col_tooltip = col_samples,
                               verbose = TRUE) {
  assertthat::assert_that(is.data.frame(.data))
  assertthat::assert_that(assertthat::is.string(col_genes))
  assertthat::assert_that(assertthat::is.string(col_samples))
  assertthat::assert_that(is.null(col_mutation_type) | assertthat::is.string(col_mutation_type))
  assertthat::assert_that(assertthat::is.string(col_tooltip))


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



  # Rank Genes based on mutation frequency / their order of appearance
  # code above already spits out genes_for_oncoplot in the appropriate order
  data_top_genes_rank <- rev(seq_along(genes_for_oncoplot))


  # Filter dataset to only include the topn/user-specified genes
  data_top_df <- .data |>
    dplyr::filter(.data[[col_genes]] %in% genes_for_oncoplot)


  # Order Genes Variable based on order
  data_top_df[[col_genes]] <- forcats::fct_relevel(data_top_df[[col_genes]], genes_for_oncoplot)

  # Sort Samples by mutated gene
  data_top_df <- data_top_df |>
    dplyr::group_by(.data[[col_samples]]) |>
    dplyr::mutate(
      SampleRankScore = score_based_on_gene_rank(mutated_genes = .data[[col_genes]], genes_informing_score = genes_for_oncoplot, gene_rank = data_top_genes_rank) # add secondary ranking based on secondary
    ) |>
    dplyr::ungroup()
  data_top_df[[col_samples]] <- forcats::fct_rev(forcats::fct_reorder(data_top_df[[col_samples]], data_top_df$SampleRankScore))

  # Consolidate to 1 row per sample-gene combo (collapse multiple mutations per gene into 1 row)
  # If col_mutation_type is supplied, will set mutation type to 'Multiple' for genes mutated multiple times in one patient
  if (!is.null(col_mutation_type)) {
    data_top_df <- data_top_df |> # Need to figure out how to fix this.
      dplyr::group_by(.data[[col_samples]], .data[[col_genes]]) |>
      dplyr::summarise(
        MutationType = unique(.data[[col_mutation_type]]) |>
          paste0(collapse = "; "),
        MutationCount = dplyr::n(),
        Tooltip = paste0(unique(.data[[col_tooltip]]), collapse = "<br>") # Edit this line to change how tooltips are collapsed
      ) |>
      dplyr::mutate(
        MutationType = ifelse(.data[["MutationCount"]] > 1, "Multi_Hit",.data[["MutationType"]])
      ) |>
      dplyr::ungroup()
  } else {
    data_top_df <- data_top_df |> # Need to figure out how to fix this.
      dplyr::group_by(.data[[col_samples]], .data[[col_genes]]) |>
      dplyr::summarise(
        MutationType = NA_character_,
        MutationCount = dplyr::n(),
        Tooltip = paste0(unique(.data[[col_tooltip]]), collapse = "<br>") # Edit this line to change how tooltips are collapsed
      ) |>
      dplyr::ungroup()
  }


  # Select just the columns we need,
  data_top_df <- data_top_df |>
   dplyr::select(
     Sample = {{ col_samples }},
     Gene = {{ col_genes }},
     MutationType = .data[["MutationType"]],
     MutationCount = .data[["MutationCount"]],
     Tooltip = .data[["Tooltip"]]
   )

  return(data_top_df)
}


# Plotting Functions ------------------------------------------------------



#' Plot oncoplot
#'
#' This function takes the output from **ggoncoplot_prep_df** and plots it.
#' Should not be exposed since it makes some assumptions about structure of input data.
#'
#' @inheritParams ggoncoplot
#' @param .data transformed data from [ggoncoplot_prep_df()] (data.frame)
#' @param margin_t,margin_r,margin_b,margin_l margin for top, right, bottom, and left side of plot. By default, unit is 'cm' but can be changed by setting `margin_unit` to any value [ggplot2::margin()] will understand (number)
#' @param margin_unit Unit of margin specification. By default is 'cm' but can be changed by setting `margin_unit` to any value [ggplot2::margin()] will understand (string)
#'
#' @inherit ggoncoplot return
#' @inherit ggoncoplot examples
ggoncoplot_plot <- function(.data,
                            show_sample_ids = FALSE,
                            palette = NULL,
                            show_ylab_title = FALSE,
                            show_xlab_title = FALSE,
                            xlab_title = "Sample",
                            ylab_title = "Gene",
                            fontsize_xlab = 16,
                            fontsize_ylab = 16,
                            fontsize_genes = 14,
                            fontsize_samples = 10,
                            tile_height = 1,
                            tile_width = 1,
                            colour_backround = "grey90",
                            colour_mutation_type_unspecified = "grey10",
                            margin_t = 0.2,
                            margin_r = 0.3,
                            margin_b = 0.2,
                            margin_l = 0.3,
                            margin_unit = "cm"
                            ) {
  check_valid_dataframe_column(.data, c("Gene", "Sample", "MutationType", "Tooltip"))

  # Invert gene factor levels
  # The gene that appears first in the levels should appear at the top of the oncoplot
  .data[["Gene"]] <- forcats::fct_rev(.data[["Gene"]])


  # Get coords of non-mutated tiles we're going to want to render in grey later
  non_mutated_tiles_df <- get_nonmutated_tiles(.data)

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
  gg <- gg +
    ggiraph::geom_tile_interactive(
      data = .data,
      ggplot2::aes_string(
        tooltip = "Tooltip",
        data_id = "Sample",
        height = tile_height,
        width = tile_width
      )
    ) +
    ggiraph::geom_tile_interactive(
      data = non_mutated_tiles_df,
      ggplot2::aes_string(
        tooltip = "Sample", # Can't just use tooltip since these don't have a value in .data. Maybe I should fix the source problem
        data_id = "Sample",
      ),
      height = tile_height,
      width = tile_width,
      fill = colour_backround
    )



  # Label axis
  gg <- gg + ggplot2::xlab(xlab_title) + ggplot2::ylab(ylab_title)

  # Add fill colour
  gg <- gg +
    ggplot2::scale_fill_manual(values = palette, na.value = colour_mutation_type_unspecified)


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

  # Show/Hide axis titles
  if(!show_xlab_title)
    gg <- gg + ggplot2::theme(axis.title.x = ggplot2::element_blank())

  if(!show_ylab_title)
    gg <- gg + ggplot2::theme(axis.title.y = ggplot2::element_blank())

  # Adjust legend position
  gg <- gg + ggplot2::theme(legend.position = "bottom")

  # Adjust legend colnumber
  gg <- gg + ggplot2::guides(fill = ggplot2::guide_legend(title = NULL, ncol = 3, keywidth=0.5))

  #Adjust legend margin
  gg <- gg + ggplot2::theme(
    legend.box.margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 5, unit = "cm")
  )

  # Adjust Margins
  gg <- gg + ggplot2::theme(
    plot.margin = ggplot2::margin(t = margin_t, r = margin_r, b = margin_b, l = margin_l, unit = margin_unit)#,
    #legend.box.margin = ggplot2::margin(t = margin_t, r = margin_r, b = margin_b, l = margin_l, unit = margin_unit)
    )

  # Adjust X scale
  gg <- gg + ggplot2::scale_x_discrete(
    drop = FALSE,
    expand = ggplot2::expansion(c(0, 0))
  )

  # Adjust Y Scale
  gg <- gg + ggplot2::scale_y_discrete(
    expand = ggplot2::expansion(c(0, 0))
  )


  return(gg)
}


# Consistent Colour Scheme
topn_to_palette <- function(.data, palette = NULL, verbose = TRUE){
  unique_impacts <- unique(.data[["MutationType"]])
  unique_impacts_minus_multiple <- unique_impacts[unique_impacts != "Multi_Hit"]

  if (all(is.na(unique_impacts))) {
    palette <- NA
  } else if (is.null(palette)) {
    mutation_dictionary <- mutationtypes::mutation_types_identify(unique_impacts_minus_multiple, verbose = verbose)

    if (mutation_dictionary == "MAF") {
      if(verbose) cli::cli_alert_success("Mutation Types are described using valid MAF terms ... using MAF palete")
        palette <- c(mutationtypes::mutation_types_maf_palette(), Multi_Hit = "black")
      palette <- palette[names(palette) %in% unique_impacts]
    } else if (mutation_dictionary == "SO") {
      if(verbose) cli::cli_alert_success("Mutation Types are described using valid SO terminology ... using SO palete")
      palette <- c(mutationtypes::mutation_types_so_palette(), Multi_Hit = "black")
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
  return(palette)
}


#' Gene barplot
#'
#' @param .data data frame output by ggoncoplot_prep_df
#' @inheritParams ggoncoplot
#' @return ggplot showing gene mutation counts
#'
#'
ggoncoplot_gene_barplot <- function(.data, fontsize_count = 14, palette = NULL, colour_mutation_type_unspecified = "grey10"){

  .data[["Gene"]] <- forcats::fct_rev(.data[["Gene"]])


  # Prepare dataframe with sample number counts
  .datacount <- dplyr::count(
      .data,
      .data[["Gene"]],
      .data[['MutationType']],
      name = "Mutations"
    ) |>
    dplyr::mutate(
      MutationType = forcats::fct_rev(
        forcats::fct_reorder(.data[["MutationType"]], .data[['Mutations']])
      )
    )

  ggplot2::ggplot(.datacount, ggplot2::aes_string(
      x = "Mutations",
      y = "Gene",
      fill = "MutationType",
      tooltip = "Mutations",
      data_id = "MutationType"
    )) +
    ggiraph::geom_col_interactive() +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      panel.grid = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm"),
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = fontsize_count)
    ) +
    ggplot2::scale_fill_manual(values = palette, na.value = colour_mutation_type_unspecified) +
    ggplot2::scale_x_continuous(position = "bottom")
}

ggoncoplot_tmb_barplot <- function(.data, col_samples, col_mutation_type, palette, colour_mutation_type_unspecified = "grey10", log10_transform = TRUE, show_ylab = FALSE,fontsize_ylab = 14, fontsize_axis_text = 11, nbreaks = 2){

  if(is.null(col_mutation_type)){
    .data[["MutationType"]] <- NA
  }
  else {
    .data <- dplyr::rename(.data, "MutationType" = {{col_mutation_type}})
  }

  df_counts <- .data |>
    dplyr::count(
      .data[[col_samples]],
      .data[["MutationType"]],
      name = "Mutations", .drop = FALSE
      )
  # Create tooltip
  df_counts$Tooltip = paste0(
    df_counts[[col_samples]], "<br>",
    "Mutations: ", df_counts[["Mutations"]]
  )



  # Main Plot
  gg <- df_counts |>
    ggplot2::ggplot(ggplot2::aes_string(y = "Mutations", x = col_samples)) +
    ggiraph::geom_col_interactive(
      ggplot2::aes(
        tooltip = .data[["Tooltip"]],
        data_id = .data[[col_samples]],
        fill = .data[["MutationType"]]
      ),
      width = 1,
      show.legend = FALSE
    )

  # Theme
  gg <- gg + ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_line(),
      axis.line.y = ggplot2::element_line(),
      axis.line.x = ggplot2::element_line(),
      panel.grid = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_text(face = "bold", size = fontsize_ylab),
      axis.text.y = ggplot2::element_text(size = fontsize_axis_text)
    )

  # Add palette arg and colour pal

  # Scales (X)
  gg <- gg + ggplot2::scale_x_discrete(drop = FALSE)

  # Scales (Y)
  trans = ifelse(log10_transform, yes = "log10", no = "identity")
  gg <- gg + ggplot2::scale_y_continuous(
    trans = trans,
    oob = scales::oob_squish_any,
    n.breaks = nbreaks,
    expand = ggplot2::expansion(c(0, 0))
    )

  # Y Axis Title
  ylabel = ifelse(log10_transform, yes = "log10\nnMuts", no = "nMuts")
  gg <- gg + ggplot2::ylab(ylabel)

  if(!show_ylab) gg <- gg + ggplot2::theme(axis.title.y = ggplot2::element_blank())

  # colour
  gg <- gg + ggplot2::scale_fill_manual(values = palette, na.value = colour_mutation_type_unspecified)
  #browser()
  return(gg)
}


#' Combine margin plots with main plot
#'
#' @param gg_main main oncoplot tileplot (ggplot)
#' @param gg_tmb barplot describing total mutations. Set to NULL to not draw barplot (ggplot)
#' @param gg_gene barplot describing number of mutated samples per gene. Set to NULL to not draw barplot (ggplot)
#'
#' @return patchwork object (or ggplot obj if both `gg_tmb` and `gg_gene` are NULL)
#'
combine_plots <- function(gg_main, gg_tmb = NULL, gg_gene = NULL){

  gg_main_margins <- gg_main$theme$plot.margin
  unit <- unique(grid::unitType(gg_main_margins))

  gg_main <- gg_main + ggplot2::theme(plot.margin = ggplot2::margin(
    t = ifelse(!is.null(gg_tmb), yes = 0, no = gg_main_margins[1]),
    r = ifelse(!is.null(gg_gene), yes = 0, no = gg_main_margins[2]),
    b = gg_main_margins[3],
    l = gg_main_margins[4],
    unit = unit
  ))

  # Both TMB and gene plots supplied
  if(!is.null(gg_tmb) & !is.null(gg_gene)){
    gg_final <- gg_tmb + patchwork::plot_spacer() + gg_main + gg_gene +
      patchwork::plot_layout(
        ncol = 2,
        widths = c(4, 1),
        heights = c(1, 10)
      )
  }
  # Only TMB
  else if(!is.null(gg_tmb) & is.null(gg_gene)){
    gg_final <- gg_tmb / gg_main +
      patchwork::plot_layout(
        heights = c(1, 10)
      )
  }
  # Only Gene
  else if(is.null(gg_tmb) & !is.null(gg_gene)){
    gg_final <- gg_main + gg_gene +
      patchwork::plot_layout(
        ncol = 2,
        widths = c(4, 1)
      )
  }
  # Neither TMB nor Gene
  else if(is.null(gg_tmb) & is.null(gg_gene)){
    gg_final <- gg_main
  }
  else
    cli::cli_abort("unexplained case when combining margin plots, package maintainer should please explicitly describe how plots should combine")


  #Add guide area down the bottom
  # gg_final <- gg_final / (patchwork::guide_area() + patchwork::plot_spacer())
  #   patchwork::plot_layout(nrow = 2, heights = c(15, 2), guides = "collect")

  return(gg_final)
}
# Utils -------------------------------------------------------------------


#' Prepare dataset for plotting
#'
#' Take a dataframe containing a column describing sample IDs (`col_sample`)
#' Filter on `col_sample` %in% samples_to_show.
#' Add any missing samples_to_show not present DF as levels of `col_sample`.
#' This way, when plotting we can use scale_x_discrete(drop=FALSE) to display all the samples we care about
#'
#'
#' @param .data dataframe with a column describing sample IDs (data.frame)
#' @param col_samples name of column in `data` containing sample IDs (character)
#' @param samples_to_show the samples we want to show in plots.
#' These samples should be the only ones represented in data.frame content,
#'  and any missing ones will be added as factor levels (character)
#'
#' @return data.frame
#'
unify_samples <- function(.data, col_samples, samples_to_show){

  # Filter to include ONLY samples in samples_to_show
  .data <- .data[.data[[col_samples]] %in% samples_to_show,]

  # Drop any extra levels based on original content
  .data[[col_samples]] <- droplevels(as.factor(.data[[col_samples]]))

  # add levels for any samples_to_show that are missing from content
  .data[[col_samples]] <- forcats::fct_expand(.data[[col_samples]], samples_to_show)

  # Ensure metadata columns are in the same order as the sequence of samples_to_show
  .data[[col_samples]] <- forcats::fct_relevel(.data[[col_samples]], samples_to_show)

  return(.data)
}

get_genes_for_oncoplot <- function(.data, col_samples, col_genes, topn, genes_to_ignore = NULL, return_extra_genes_if_tied = FALSE, genes_to_include = NULL, verbose = TRUE){
  # Look exclusively at a custom set of genes
  if (!is.null(genes_to_include)) {
    genes_not_found <- genes_to_include[!genes_to_include %in% .data[[col_genes]]]

    if (length(genes_not_found) == length(genes_to_include)) {
      cli::cli_abort("Couldn't find any of the genes you supplied in your dataset. Either no samples have mutations in these genes, or you've got the wrong gene names")
    }

    if (length(genes_not_found) > 0) {
      if(verbose){
        cli::cli_warn(
          c(
            "Failed to find the following [{length(genes_not_found)}] genes in your dataset",
            ">" = "{genes_not_found}",
            "!" = "Either no samples have mutations in the above genes, or you've got the wrong gene names"
            )
        )
        #cli::cli_alert("{genes_not_found}")
        #cli::cli_alert_warning("Either no samples have mutations in the above genes, or you've got the wrong gene names")
      }
      # Filter out genes that aren't found
      genes_to_include <- genes_to_include[!genes_to_include %in% genes_not_found]
    }

    # filter out any 'genes_to_ignore'
    genes_to_include <- genes_to_include[!genes_to_include %in% genes_to_ignore]
    genes_for_oncoplot <- genes_to_include
  }
  # Look only at the topn mutated genes
  else{
    genes_for_oncoplot <- identify_topn_genes(
      .data = .data,
      col_samples = col_samples,
      col_genes = col_genes,
      topn = topn,
      return_extra_genes_if_tied = return_extra_genes_if_tied,
      genes_to_ignore = genes_to_ignore,
      verbose = verbose
    )
  }
}

#' Identify top genes from a mutation df
#'
#' Identify top genes from a mutation df
#'
#' @inheritParams ggoncoplot
#'
#' @return vector of topn genes. Their order will be their rank (most mutated = first) (character)
#'
identify_topn_genes <- function(.data, col_samples, col_genes, topn, genes_to_ignore = NULL, return_extra_genes_if_tied = FALSE, verbose = TRUE){
  assertthat::assert_that(assertthat::is.flag(return_extra_genes_if_tied))
  assertthat::assert_that(is.null(genes_to_ignore) | is.character(genes_to_ignore))
  assertthat::assert_that(assertthat::is.number(topn))
  assertthat::assert_that(topn > 0)
  assertthat::assert_that(assertthat::is.flag(verbose))

  # Identify top genes by frequency
  df_data_gene_counts <- .data |>
    dplyr::ungroup() |>
    dplyr::distinct(.data[[col_samples]], .data[[col_genes]]) |> # This line stops multiple mutations of the same gene in the same sample counting multiple times towards the mutation frequency.
    dplyr::count(.data[[col_genes]]) |>
    dplyr::arrange(dplyr::desc(.data[["n"]]), .data[[col_genes]]) # Sort by count, then by alphabetic order (So if return_extra_genes_if_tied == FALSE we break ties based on alphabetic order, not order of appearence in DF)

  if(!is.null(genes_to_ignore)){
    df_data_gene_counts <- df_data_gene_counts |>
      dplyr::filter(! .data[[col_genes]] %in% genes_to_ignore)
  }

  # Would have to add |> collect() somewhere in this function to support DBs

  data_top_genes_df <- df_data_gene_counts |>
    dplyr::slice_max(.data$n, n = topn, with_ties = return_extra_genes_if_tied) # Set with_ties = TRUE to allow topn genes to return extra genes if there are ties in # of samples mutated

  top_genes <- data_top_genes_df[[col_genes]]


  if(length(top_genes) > topn){
    if(verbose) cli::cli_alert_info(
      'The top [{topn}] genes were requested.
      BUT due to ties, the top {length(top_genes)} genes were returned.
      Change this behaviour by setting {.arg return_extra_genes_if_tied = FALSE}'
      )
  }
  else if(length(top_genes) < topn){
    if(verbose) cli::cli_alert_info(
      'The top [{topn}] genes were requested.
      BUT not enough genes are in data, so all [{length(top_genes)}] genes in dataset were returned'
    )
  }

  return(top_genes)
}


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
      axis.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::unit(c(0, 0, 0, 0), "cm")
    )
}


#' data.frame has colnames
#'
#' Assert that data.frame contains a set of user defined column names.
#'
#' data.frame may have any additional colnames.
#' It just has to have AT LEAST the columns specified by `colnames`
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


#' Get data.frame o
#'
#' Takes same .data input as ggoncoplot and returns a dataframe with 'Sample' and 'Gene' columns
#' ONLY for sample-gene pairs that are unmutated. This lets us colour render them separately (as grey)
#'
#' @inheritParams ggoncoplot_plot
#'
#' @return  a dataframe with 'Sample' and 'Gene' columns ONLY for sample-gene pairs that are unmutated. This lets us colour render them separately (as grey)  (data.frame)
get_nonmutated_tiles <- function(.data){
  samples <- levels(.data[['Sample']])

  non_mutated_tiles_df  <- expand.grid(
    Sample = samples,
    Gene = unique(.data[["Gene"]])
  )

  nomutations <- ! paste0(
    non_mutated_tiles_df[['Sample']],
    non_mutated_tiles_df[['Gene']]
  ) %in%
    unique(
      paste0(
        .data[['Sample']],
        .data[["Gene"]]
      )
    )

  non_mutated_tiles_df[nomutations,]
}
