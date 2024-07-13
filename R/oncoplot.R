
# Globals -----------------------------------------------------------------
utils::globalVariables(
  c(
    "Gene", "MutationType", "Pathway", "Sample", "Tooltip", "MutationCount",
    "Mutations", "count", "fill", "y", ".data"
  )
)

# Oncoplot ----------------------------------------------------------------

#' ggoncoplot
#' @importFrom patchwork plot_layout
#'
#' @param data data for oncoplot. A data.frame with 1 row per mutation in your cohort. Must contain columns describing gene_symbols and sample_identifiers (data.frame)
#' @param col_genes name of **data** column containing gene names/symbols (string)
#' @param col_samples name of **data** column containing sample identifiers (string)
#' @param col_mutation_type name of **data** column describing mutation types (string, optional)
#' @param genes_to_include specific genes to include in the oncoplot (character, optional)
#' @param genes_to_ignore names of the genes that should be ignored (character, optional)
#' @param col_tooltip name of **data** column containing whatever information you want to display in (string, defaults to col_samples)
#' @param topn how many of the top genes to visualize. Ignored if `genes_to_include` is supplied (number, default 10)
#' @param return_extra_genes_if_tied instead of strictly returning `topn` genes,
#' in the case of ties (where multiple genes are mutated in the exact same number of samples, complicating selection of top n genes), return all tied genes (potentially more than topn).
#' If FALSE, will return strictly `topn` genes, breaking ties based on order of appearance in dataset (flag, default FALSE)
#' @param draw_gene_barplot add a barplot describing number of samples with each gene mutated (right side) (flag, default FALSE)
#' @param draw_tmb_barplot add a barplot describing total number of mutations in each sample (above main plot). If a single gene is mutated multiple times, all mutations are counted towards total (flag, default FALSE)
#' @param copy value to copy to clipboard when an oncoplot tile is clicked (string, one of 'sample', 'gene', 'tooltip', 'mutation_type', 'nothing', default 'sample')
#' @param palette a named vector mapping all possible mutation types (vector names) to colors (vector values, optional)
#' @param metadata dataframe describing sample level metadata.
#' One column must contain unique sample identifiers. Other columns can describe numeric / categorical metadata (data.frame, optional)
#' @param metadata_palette A list of named vectors. List names correspond to metadata column names (categorical only). Vector names to levels of columns. Vector values are colors, the vector names are used to map values in data to a color. (optional)
#' @param col_samples_metadata which column in metadata data.frame describes sample identifiers (string, defaults to col_samples)
#' @param cols_to_plot_metadata names of columns in metadata that should be plotted (character, optional)
#' @param metadata_require_mutations filter out samples from metadata lacking any mutations in data (flag, default TRUE)
#' @param pathway a two column dataframe describing pathway. The column containing gene names should have the same name as \strong{col_gene} (data.frame, optional)
#' @param col_genes_pathway which column in pathway data.frame describes gene names (string, defaults to col_genes)
#' @param show_all_samples show all samples in oncoplot, even if they don't have mutations in the selected genes. Samples only described in metadata but with no mutations at all are still filtered out by default, but you can show these too by setting `metadata_require_mutations = FALSE` (flag, default FALSE)
#' @param total_samples Strategy for calculating the total number of samples.
#' This value is used to compute the proportion of mutation recurrence displayed in the tooltip when hovering over the gene barplot,
#' or as a text annotation when \code{ggoncoplot_options(show_genebar_labels = TRUE)} is set to TRUE.
#'
#' Possible values:
#' \itemize{
#'   \item \strong{any_mutations}: All the samples that are in \code{data} (the mutation dataset), irrespective of whether they are on the oncoplot or not.
#'   \item \strong{oncoplot}: Only the samples that are present on the oncoplot.
#'   \item \strong{all}: All the samples in either \code{data} or \code{metadata}.
#' }
#' @param interactive should plot be interactive (boolean, default TRUE)
#' @param verbose verbose mode (flag, default TRUE)
#' @param options a list of additional visual parameters created by calling [ggoncoplot_options()]. See \code{\link{ggoncoplot_options}} for details.
#'
#' @inheritDotParams gg1d::gg1d
#'
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
#' gbm_clinical_csv <- system.file(
#'   package = "ggoncoplot",
#'   "testdata/GBM_tcgamutations_mc3_clinical.csv"
#' )
#'
#' gbm_df <- read.csv(file = gbm_csv, header = TRUE)
#' gbm_clinical_df <- read.csv(file = gbm_clinical_csv, header = TRUE)
#'
#' # Plot Basic Oncoplot
#' ggoncoplot(
#'   gbm_df,
#'   "Hugo_Symbol",
#'   "Tumor_Sample_Barcode",
#'   col_mutation_type = "Variant_Classification",
#'   metadata = gbm_clinical_df,
#'   cols_to_plot_metadata = "gender"
#' )
#'
#' # Customise how the Oncoplot looks
#' ggoncoplot(
#'   gbm_df,
#'   "Hugo_Symbol",
#'   "Tumor_Sample_Barcode",
#'   col_mutation_type = "Variant_Classification",
#'   metadata = gbm_clinical_df,
#'   cols_to_plot_metadata = "gender",
#'
#'   # Customise Visual Options
#'   options = ggoncoplot_options(
#'       xlab_title = "Glioblastoma Samples",
#'       ylab_title = "Top 10 mutated genes"
#'   )
#' )
ggoncoplot <- function(data,
                       col_genes, col_samples,
                       col_mutation_type = NULL,
                       genes_to_include = NULL,
                       genes_to_ignore = NULL,
                       col_tooltip = col_samples,
                       topn = 10,
                       return_extra_genes_if_tied = FALSE,
                       draw_gene_barplot = FALSE,
                       draw_tmb_barplot = FALSE,
                       copy = c('sample', 'gene', 'tooltip', 'mutation_type', 'nothing'),
                       palette = NULL,
                       metadata = NULL,
                       metadata_palette = NULL,
                       col_samples_metadata = col_samples,
                       cols_to_plot_metadata = NULL,
                       metadata_require_mutations = TRUE,
                       pathway = NULL,
                       col_genes_pathway = col_genes,
                       show_all_samples = FALSE,
                       total_samples = c('any_mutations', 'all', 'oncoplot'),
                       interactive = TRUE,
                       options = ggoncoplot_options(),
                       verbose = TRUE,
                       ...
                       ) {


  # Assertions --------------------------------------------------------------
  assertions::assert_dataframe(data)
  assertions::assert(nrow(data) > 0)
  assertions::assert_string(col_genes)
  assertions::assert_string(col_samples)
  if(!is.null(genes_to_include)) assertions::assert_character(genes_to_include)
  assertions::assert_string(col_tooltip)
  assertions::assert_number(topn)
  assertions::assert_flag(verbose)
  assertions::assert_flag(draw_gene_barplot)
  assertions::assert_flag(draw_tmb_barplot)


  if(!is.null(metadata)){
    assertions::assert_dataframe(metadata)
    assertions::assert_names_include(metadata, col_samples_metadata)
    assertions::assert_no_duplicates(metadata[[col_samples_metadata]], arg_name = "Metadata Sample Column")
  }
  if(!is.null(col_mutation_type)){
    assertions::assert_string(col_mutation_type)
    # Assert mutation type is a valid column name
    assertions::assert_names_include(data, col_mutation_type)

    # If column type is a factor, convert to character
    if(is.factor(data[[col_mutation_type]])) data[[col_mutation_type]] <- as.character(data[[col_mutation_type]])

    # Assert column type = character
    #assertions::assert_character(data[[col_mutation_type]])
    assertions::assert_no_missing(data[[col_mutation_type]], arg_name = paste0("Mutation Type Column: ", col_mutation_type))
    assertions::assert_excludes(data[[col_mutation_type]], illegal = "", msg = "{.strong Mutation Type} column cannot contain zero-length strings")
  }

  if(!is.null(pathway)){
    assertions::assert_dataframe(pathway)
    assertions::assert_string(col_genes_pathway)
    assertions::assert_names_include(pathway, col_genes_pathway)
    assertions::assert_character(pathway[[col_genes_pathway]])


    assertions::assert(ncol(pathway) == 2, msg = "Pathway dataframe must have exactly 2 columns")
    col_pathways_pathway <- (colnames(pathway)[!colnames(pathway) %in% col_genes_pathway])[1]
    if(verbose) message("Found pathway column: ", col_pathways_pathway)
    assertions::assert_character(pathway[[col_pathways_pathway]])
    assertions::assert_excludes(pathway[[col_pathways_pathway]], illegal = "Other", msg = "You have a pathway named 'Other' in your pathway data. This is not allowed because ggoncoplot automaticaly sets all genes without a pathway, to pathway 'Other'. To fix, simply remove all rows where pathway = 'Other'")

    assertions::assert_no_missing(pathway[[col_genes_pathway]])
    assertions::assert_no_duplicates(pathway[[col_genes_pathway]])
    assertions::assert_no_missing(pathway[[col_pathways_pathway]])

    # Reorder columns so pathway[[1]] gives you genes and pathway[[2]] gives you pathways
    pathway <- pathway[c(col_genes_pathway, col_pathways_pathway)]
  }

  #Assert sample column sensible
  if(is.factor(data[[col_samples]])) data[[col_samples]] <- as.character(data[[col_samples]])
  assertions::assert_character(data[[col_samples]])
  assertions::assert_no_missing(data[[col_samples]])
  assertions::assert_excludes(data[[col_samples]], illegal = "", msg = "{.strong Sample} column cannot contain zero-length strings") # Asserts no empty string

  #Assert gene column sensible
  if(is.factor(data[[col_genes]])) data[[col_genes]] <- as.character(data[[col_genes]])
  assertions::assert_character(data[[col_genes]])
  assertions::assert_no_missing(data[[col_genes]])
  assertions::assert_excludes(data[[col_genes]], illegal = "", msg = "{.strong Gene} column cannot contain zero-length strings") # Asserts no empty string

  # Assert options are produced by ggoncoplot_options()
  assertions::assert_class(options, "ggoncoplot_options")

  # Argument matching
  copy <- rlang::arg_match(copy)
  total_samples <- rlang::arg_match(total_samples)

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
    lgl_samples_have_muts <- metadata[[col_samples_metadata]] %in% unique(data[[col_samples]])
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
    data = data,
    col_samples = col_samples,
    col_genes = col_genes,
    topn = topn,
    genes_to_ignore = genes_to_ignore,
    return_extra_genes_if_tied = return_extra_genes_if_tied,
    genes_to_include = genes_to_include,
    verbose = verbose
  )

  # Rerank genes based on pathway data.frame
  if(!is.null(pathway)){
    genes_for_oncoplot <- rank_genes_based_on_pathways(
      gene_pathway_map = pathway,
      generanks = genes_for_oncoplot,
      pathwayranks = unique(pathway[[2]])
    )
  }

  # Preprocess dataframe ----------------------------------------------------
  # Get dataframe with 1 row per sample-gene pair
  #TODO: add pathway argument, and ensure pathway variable gets returned (factor with levels reasonably sorted)
  data_top_df <- ggoncoplot_prep_df( # Add a samples_for_oncoplot
    data = data,
    col_genes = col_genes, col_samples = col_samples,
    col_mutation_type = col_mutation_type,
    col_tooltip = col_tooltip,
    pathway = pathway,
    genes_for_oncoplot = genes_for_oncoplot,
    verbose=verbose
  )


  # Sample order ----------------------------------------------
  # Get Sample Order,
  samples_with_mutations_in_selected_genes <- levels(droplevels(data_top_df[["Sample"]]))
  samples_with_any_mutations <- unique(data[[col_samples]])

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
  data <- unify_samples(data = data, col_samples = col_samples, samples_to_show = samples_to_show)
  data_top_df <- unify_samples(data = data_top_df, col_samples = "Sample", samples_to_show = samples_to_show)
  metadata <- unify_samples(data = metadata, col_samples = col_samples_metadata, samples_to_show = samples_to_show)

  # Samples in Onocplot
  samples_in_oncoplot = data_top_df$Sample

  # Calculate the total number of samples (useful for calculations of prevl
  n_total_samples <-
    if(total_samples == "any_mutations")
      n_total_samples <- dplyr::n_distinct(samples_with_any_mutations)
    else if(total_samples == "all")
      dplyr::n_distinct(all_sample_ids)
    else if(total_samples == "oncoplot")
      dplyr::n_distinct(samples_in_oncoplot)
    else
      stop("the ggoncoplot package developer messed up. Unnacounted value of total_samples: ", total_samples)

  # Palette -----------------------------------------------------------------
  palette <- topn_to_palette(data = data_top_df, palette = palette, verbose = verbose)



  # Draw main plot --------------------------------------------------------
  gg_main <- ggoncoplot_plot(
    data = data_top_df,
    show_sample_ids = options$show_sample_ids,
    palette = palette,
    xlab_title = options$xlab_title,
    ylab_title = options$ylab_title,
    fontsize_xlab = options$fontsize_xlab,
    fontsize_ylab = options$fontsize_ylab,
    fontsize_genes = options$fontsize_genes,
    fontsize_samples = options$fontsize_samples,
    fontsize_legend_text = options$fontsize_legend_text,
    fontsize_legend_title = options$fontsize_legend_title,
    legend_key_size = options$legend_key_size,
    copy = copy,
    tile_height = options$tile_height,
    tile_width = options$tile_width,
    colour_backround = options$colour_backround,
    margin_t = margin_main_t,
    margin_r = margin_main_r,
    margin_b = margin_main_b,
    margin_l = margin_main_l,
    legend_title = if(is.null(col_mutation_type)) "Mutation Type" else beautify(col_mutation_type),
    show_ylab_title = options$show_ylab_title,
    show_xlab_title = options$show_ylab_title,
    margin_unit = margin_units,
    colour_mutation_type_unspecified = options$colour_mutation_type_unspecified,
    colour_pathway_text = options$colour_pathway_text,
    colour_pathway_bg = options$colour_pathway_bg,
    colour_pathway_outline = options$colour_pathway_outline,
    pathway_text_angle = options$pathway_text_angle,
    fontsize_pathway = options$fontsize_pathway,
    ggoncoplot_guide_ncol = options$ggoncoplot_guide_ncol,
    show_legend_titles = options$show_legend_titles
  )

  # Draw marginal plots -----------------------------------------------------
  gg_gene_barplot = NULL
  gg_tmb_barplot = NULL
  gg_metadata = NULL

  ## Gene Barplot -----------------------------------------------------------
  if(draw_gene_barplot){
    gg_gene_barplot <- ggoncoplot_gene_barplot(
      data = data_top_df,
      total_samples = n_total_samples,
      palette = palette,
      fontsize_count = options$fontsize_count,
      colour_mutation_type_unspecified = options$colour_mutation_type_unspecified,
      show_axis = options$show_axis_gene,
      genebar_label_padding = options$genebar_label_padding,
      show_genebar_labels = options$show_genebar_labels,
      genebar_label_nudge = options$genebar_label_nudge,
      only_pad_if_labels_shown = options$genebar_only_pad_when_labels_shown,
      digits_to_round_to = options$genebar_label_round,
      genebar_scale_breaks = options$genebar_scale_breaks,
      genebar_scale_n_breaks = options$genebar_scale_n_breaks
    )

  }

  ## TMB plot  -----------------------------------------------------------
  if(draw_tmb_barplot){
    gg_tmb_barplot <- ggoncoplot_tmb_barplot(
      data = data,
      col_samples = col_samples,
      col_mutation_type = col_mutation_type,
      log10_transform = options$log10_transform_tmb,
      fontsize_ylab = options$fontsize_tmb_title,
      fontsize_axis_text = options$fontsize_tmb_axis,
      show_ylab = options$show_ylab_title_tmb,
      palette = palette,
      colour_mutation_type_unspecified = options$colour_mutation_type_unspecified,
      scientific = options$scientific_tmb,
      show_axis = options$show_axis_tmb,
      verbose = verbose
    )

  }

  ## Draw sample metadata plots ---------------------------------------------------------
  if(!is.null(metadata)){
    gg_metadata <- gg1d::gg1d(
      metadata,
      col_id = col_samples_metadata,
      cols_to_plot = cols_to_plot_metadata,
      interactive = FALSE,
      verbose = if(verbose) 1 else 0,
      cli_header = "Plotting Sample Metadata",

      # Tile width
      width = options$tile_width,

      # Axis Text Fontsizes
      fontsize_y_text = options$fontsize_metadata_text,

      # Legend Fontsizes
      legend_title_size = options$fontsize_metadata_legend_title,
      legend_text_size = options$fontsize_metadata_legend_text,
      fontsize_barplot_y_numbers = options$fontsize_metadata_barplot_y_numbers,

      # Legend Layout
      show_legend_titles = options$show_legend_titles, #default TRUE
      legend_nrow = options$metadata_legend_nrow, # default NULL
      legend_ncol = options$metadata_legend_ncol, # default NULL
      legend_key_size = options$metadata_legend_key_size,

      # Dealing with NAs
      na_marker = options$metadata_na_marker,
      na_marker_size = options$metadata_na_marker_size,

      # Processing Levels
      maxlevels = options$metadata_maxlevels,

      # Numeric Values
      numeric_plot_type = options$metadata_numeric_plot_type,
      legend_orientation_heatmap = options$metadata_legend_orientation_heatmap,

      # Palettes
      palettes = metadata_palette,
      # colours_default = metadata_colours_default,
      # colours_default_logical = metadata_colours_default_logical,
      # colours_missing = metadata_colours_missing,
      y_axis_position = "left"
    )
  }

  ## Combine marginal plots -----------------------------------------------------------
  gg_final <- combine_plots(
    gg_main,
    gg_tmb = gg_tmb_barplot,
    gg_gene = gg_gene_barplot,
    gg_metadata = gg_metadata,
    gg_tmb_height = options$plotsize_tmb_rel_height,
    gg_gene_width = options$plotsize_gene_rel_width,
    gg_metadata_height = options$plotsize_metadata_rel_height
    )


  ## Control Look of oncoplot + marginal plots
  if(!options$show_legend){
   gg_final <- gg_final & ggplot2::theme(legend.position = "none")
  }


  # Make Interactive -------------------------------------------------------
  # Turn gg into an interactive ggiraph object if interactive = TRUE
  if (interactive) {
    gg_final <- ggiraph::girafe(
      width_svg = options$interactive_svg_width, height_svg = options$interactive_svg_height,
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
#' @param data data for oncoplot. A data.frame with 1 row per mutation in your cohort. Must contain columns describing gene_symbols and sample_identifiers (data.frame)
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
#'   data = gbm_df,
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
ggoncoplot_prep_df <- function(data,
                               col_genes,
                               col_samples,
                               genes_for_oncoplot,
                               col_mutation_type = NULL,
                               col_tooltip = col_samples,
                               pathway = NULL,
                               verbose = TRUE) {
  assertions::assert_dataframe(data)
  assertions::assert_string(col_genes)
  assertions::assert_string(col_samples)
  if(!is.null(col_mutation_type)) assertions::assert_string(col_mutation_type)
  assertions::assert_string(col_tooltip)


  # Check specified columns are in data
  data_colnames <- names(data)

  check_valid_dataframe_column(
    data = data,
    colnames = c(
      col_samples,
      col_genes,
      col_tooltip
    )
  )

  # Check optional columns are in data
  if (!is.null(col_mutation_type)) {
    check_valid_dataframe_column(data = data, colnames = col_mutation_type)
  }

  # Ensure Sample Column is A factor
  data[[col_samples]] <- as.factor(data[[col_samples]])



  # Rank Genes based on mutation frequency / their order of appearance
  # code above already spits out genes_for_oncoplot in the appropriate order
  data_top_genes_rank <- rev(seq_along(genes_for_oncoplot))

  # Rank genes based on pathway

  # Filter dataset to only include the topn/user-specified genes
  data_top_df <- data |>
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

  # Add sample Identifier to top of tooltip if not already the tooltip
  if(col_tooltip != col_samples)
    data_top_df[["Tooltip"]] <- paste0("<strong>", data_top_df[[col_samples]], "</strong><br>", data_top_df[["Tooltip"]])


  # Select just the columns we need,
  data_top_df <- data_top_df |>
   dplyr::select(
     Sample = {{ col_samples }},
     Gene = {{ col_genes }},
     MutationType = MutationType,
     MutationCount = MutationCount,
     Tooltip = Tooltip
   )

  # Add pathway column
  if(!is.null(pathway)){
    # Create Pathway Column
    data_top_df[["Pathway"]] <- pathway[[2]][match(data_top_df[["Gene"]], pathway[[1]])]
    data_top_df[["Pathway"]] <- ifelse(is.na(data_top_df[["Pathway"]]), "Other", data_top_df[["Pathway"]])
    data_top_df[["Pathway"]] <- as.factor(data_top_df[["Pathway"]])

    # Sort based on order of appearance in pathway df
    data_top_df[["Pathway"]] <- forcats::fct_relevel(.f = data_top_df[["Pathway"]], unique(pathway[[2]]))

    # TODO: Add an alternate sort based on samples mutated and an argument that lets the user choose sorting style

  }

  return(data_top_df)
}


# Plotting Functions ------------------------------------------------------



#' Plot oncoplot
#'
#' This function takes the output from **ggoncoplot_prep_df** and plots it.
#' Should not be exposed since it makes some assumptions about structure of input data.
#'
#' @inheritParams ggoncoplot
#' @inheritParams ggoncoplot_options
#' @param data transformed data from [ggoncoplot_prep_df()] (data.frame)
#' @param margin_t,margin_r,margin_b,margin_l margin for top, right, bottom, and left side of plot. By default, unit is 'cm' but can be changed by setting `margin_unit` to any value [ggplot2::margin()] will understand (number)
#' @param margin_unit Unit of margin specification. By default is 'cm' but can be changed by setting `margin_unit` to any value [ggplot2::margin()] will understand (string)
#' @param legend_title name of legend title (string)
#' @inherit ggoncoplot return
#' @inherit ggoncoplot examples
ggoncoplot_plot <- function(data,
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
                            fontsize_legend_title = 12,
                            fontsize_legend_text = 12,
                            tile_height = 1,
                            tile_width = 1,
                            copy = c('sample', 'gene', 'tooltip', 'mutation_type', 'nothing'),
                            colour_backround = "grey90",
                            colour_mutation_type_unspecified = "grey10",
                            fontsize_pathway = 16,
                            colour_pathway_text = "white",
                            colour_pathway_bg = "grey10",
                            colour_pathway_outline = "black",
                            pathway_text_angle = 0,
                            legend_title = "Mutation Type",
                            show_legend_titles = TRUE,
                            ggoncoplot_guide_ncol = 2,
                            legend_key_size = 0.3,
                            margin_t = 0.2,
                            margin_r = 0.3,
                            margin_b = 0.2,
                            margin_l = 0.3,
                            margin_unit = "cm"
                            ) {
  copy <- rlang::arg_match(copy)
  check_valid_dataframe_column(data, c("Gene", "Sample", "MutationType", "Tooltip"))

  # Invert gene factor levels
  # The gene that appears first in the levels should appear at the top of the oncoplot
  data[["Gene"]] <- forcats::fct_rev(data[["Gene"]])


  # Get coords of non-mutated tiles we're going to want to render in grey later
  non_mutated_tiles_df <- get_nonmutated_tiles(data)

  # Figure out which colum name to copy on click
  copy_column <- dplyr::case_match(
    copy,
    "sample" ~ "Sample",
    "gene" ~ "Gene",
    "tooltip" ~ "Tooltip",
    "mutation_type" ~ "MutationType",
    .default = NA_character_
    )

  # Create ggplot
  gg <- ggplot2::ggplot(
    data = data,
    mapping = ggplot2::aes(
      y = Gene,
      x = Sample,
      fill = MutationType
    )
  )

  # Add interactive/non-interactive geom layer
  gg <- gg +
    ggiraph::geom_tile_interactive(
      data = data,
      ggplot2::aes(
        tooltip = Tooltip,
        data_id = Sample,
        onclick = if (copy == "nothing") NULL else paste0('navigator.clipboard.writeText("',.data[[copy_column]],'")'),
        height = {{tile_height}},
        width = {{tile_width}}
      )
    ) +
    ggiraph::geom_tile_interactive(
      data = non_mutated_tiles_df,
      ggplot2::aes(
        tooltip = Sample, # Can't just use tooltip since these don't have a value in .data. Maybe I should fix the source problem
        data_id = Sample,
      ),
      height = tile_height,
      width = tile_width,
      fill = colour_backround
    )

  # Facet by pathway
  if("Pathway" %in% colnames(data)){
    gg <- gg + ggiraph::facet_grid_interactive(
      rows = ggplot2::vars(Pathway),
      scales = "free_y",
      space = "free_y",
      switch = "y",
      labeller = ggiraph::labeller_interactive(ggplot2::aes(tooltip = Pathway)),
      interactive_on = "both"
    )
  }

  # Label axis
  gg <- gg + ggplot2::xlab(xlab_title) + ggplot2::ylab(ylab_title)

  # Add fill colour
  gg <- gg +
    ggplot2::scale_fill_manual(values = palette, na.value = colour_mutation_type_unspecified)


  # Apply default theme
  gg <- gg + theme_oncoplot_default(
    show_legend_titles = show_legend_titles,
    fontsize_legend_title = fontsize_legend_title,
    fontsize_legend_text = fontsize_legend_text
    )

  # Add line between genes
  gg <- gg + ggplot2::geom_hline(yintercept = seq(0, length(unique(data[['Gene']]))) + .5, color="gray30")

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
  gg <- gg + ggplot2::theme(legend.position = "right")

  # Adjust legend colnumber (and set title)
  if(!show_legend_titles) legend_title <- NULL
  gg <- gg + ggplot2::guides(fill = ggplot2::guide_legend(title = legend_title, ncol = ggoncoplot_guide_ncol, keywidth=legend_key_size, title.hjust = 0))

  #Adjust legend margin
  gg <- gg + ggplot2::theme(
    legend.box.margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 5, unit = "cm")
  )

  # Adjust Margins
  gg <- gg + ggplot2::theme(
    plot.margin = ggplot2::margin(t = margin_t, r = margin_r, b = margin_b, l = margin_l, unit = margin_unit)
  )

  # Adjust X scale
  gg <- gg + ggplot2::scale_x_discrete(
    drop = FALSE,
    expand = ggplot2::expansion(c(0, 0))
  )

  # Adjust Y Scale
  gg <- gg + ggplot2::scale_y_discrete(
    expand = ggplot2::expansion(c(0, 0)), position = "left"
  )

  # Adjust pathway facet properties
  gg <- gg + ggplot2::theme(
    strip.text.y.left =  ggiraph::element_text_interactive(
      size = fontsize_pathway,
      angle = pathway_text_angle,
      color = colour_pathway_text,
      face = "bold"
      ),

    #strip.text.y = ggplot2::element_blank(),
    strip.placement = "outside",strip.clip = "on",
    strip.background = ggplot2::element_rect(fill = colour_pathway_bg, colour = colour_pathway_outline)
    #strip.background = ggplot2::element_blank()
    )

  return(gg)
}


# Consistent Colour Scheme
topn_to_palette <- function(data, palette = NULL, verbose = TRUE){
  unique_impacts <- unique(data[["MutationType"]])
  unique_impacts_minus_multiple <- unique_impacts[unique_impacts != "Multi_Hit"]
  #browser()
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
      if(any(grepl(pattern = "&", x = unique_impacts, fixed = TRUE))) cli::cli_abort("Found ampersand (&) delimited SO mutation impacts. Please run {.code mutationtypes::select_most_severe_consequence_so()} on your mutation_type column before feeding data into ggoncoplot")
    } else { # What if hits don't map well to
      cli::cli_h1("Variant Type Ontology Unknown")
      cli::cli_alert_warning("Mutation Types are not perfectly described with any known ontology.
                               Using an RColorBrewer palette by default.
                               When running this plot with other datasets, it is possible the colour scheme may differ.
                               We {.strong STRONGLY reccomend} supplying a custom MutationType -> colour mapping using the {.arg palette} argument")

      # data[['MutationType']] <- forcats::fct_infreq(f = data[['MutationType']])
      rlang::check_installed("RColorBrewer", reason = "To create default palette for `ggoncoplot()`")
      if(length(unique_impacts) > 12){
        cli::cli_abort("Too many unique Mutation Types for automatic palette generation (need <=12, not {length(unique_impacts)}). Please supply a custom Mutation Type -> colour mapping using the {.arg palette} argument")
      }
      palette <- RColorBrewer::brewer.pal(n = 12, name = "Paired")
    }
  } else { # What if custom palette is supplied?
    if (!all(unique_impacts %in% names(palette))) {
      terms_without_mapping <- unique_impacts[!unique_impacts %in% names(palette)]
      cli::cli_abort("Please add colour mappings for the following terms: {terms_without_mapping}")
      palette <- palette[names(palette) %in% unique_impacts]
    }
  }
  #browser()
  return(palette)
}


#' Gene barplot
#'
#' @param data data frame output by ggoncoplot_prep_df
#' @param show_axis show axis text/ticks/line (flag)
#' @param only_pad_if_labels_shown should expansion to x axis be applied if bar labels aren't shown?
#' @param digits_to_round_to how many digits to round recurrence proportions to
#' @inheritParams ggoncoplot
#' @inheritParams ggoncoplot_options
#' @return ggplot showing gene mutation counts
#'
ggoncoplot_gene_barplot <- function(data, fontsize_count = 14, palette = NULL,
                                    colour_mutation_type_unspecified = "grey10",
                                    show_axis, total_samples,
                                    show_genebar_labels = TRUE,
                                    genebar_label_nudge = 2,
                                    genebar_label_padding = 0.2,
                                    only_pad_if_labels_shown = TRUE,
                                    digits_to_round_to = 0,
                                    genebar_scale_n_breaks = 3,
                                    genebar_scale_breaks = ggplot2::waiver()
                                    ){

  data[["Gene"]] <- forcats::fct_rev(data[["Gene"]])

  if(!show_genebar_labels & only_pad_if_labels_shown) genebar_label_padding <- 0

  # Main plot
  gg <- ggplot2::ggplot(data, ggplot2::aes(
      y = Gene,
    )) +
    ggiraph::geom_bar_interactive(
      ggplot2::aes(
        fill = MutationType,
        data_id = Gene,
        tooltip = paste0(
          "Total Samples Mutated: ", ggplot2::after_stat(stats::ave(count, y, FUN = sum)),
          " (", as_pct(ggplot2::after_stat(stats::ave(count, y, FUN = sum)/total_samples), digits = digits_to_round_to) ," of all samples)",
          "<br/>",
          ggplot2::after_stat(fill),": ", ggplot2::after_stat(count),
          " (", as_pct(ggplot2::after_stat(count / stats::ave(count, y, FUN = sum)), digits = digits_to_round_to), " of all mutations in this gene)"
        ),
        ),
      stat="count"
      ) +
    ggplot2::coord_cartesian(clip="off") +
    ggplot2::scale_y_discrete(expand = ggplot2::expansion(c(0, 0))) +
    ggplot2::scale_x_continuous(
      position = "bottom",
      expand = ggplot2::expansion(mult = c(0, genebar_label_padding)),
      breaks = genebar_scale_breaks,
      n.breaks = genebar_scale_n_breaks,
    )

  if(show_genebar_labels)
    gg <- gg +  ggplot2::geom_text(
      ggplot2::aes(label = as_pct(ggplot2::after_stat(stats::ave(count, y, FUN = sum) / total_samples),digits = digits_to_round_to),  group = Gene),
      stat = 'count',
      #fill = "white",label.size = 0, alpha = 0.5,
      hjust=0,
      nudge_x = genebar_label_nudge
    )

  # Facet by Pathway
  if(!is.null(data[["Pathway"]])){
   gg <- gg + ggplot2::facet_grid(
     rows = ggplot2::vars(Pathway),
     scales = "free_y", space = "free_y",
   )
  }

  # Theming
    gg <- gg +
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
        axis.text.x = ggplot2::element_text(size = fontsize_count),
        strip.text.y = ggplot2::element_blank(),
        strip.background = ggplot2::element_blank()
      )

  # Add colours
  gg <- gg +
    ggplot2::scale_fill_manual(values = palette, na.value = colour_mutation_type_unspecified)

  # Show / hide main axis
  if(!show_axis){
    gg <- gg + ggplot2::theme(
      axis.line.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
      )
  }

  return(gg)
}

ggoncoplot_tmb_barplot <- function(data, col_samples, col_mutation_type, palette, colour_mutation_type_unspecified = "grey10", log10_transform = TRUE, show_ylab = FALSE,fontsize_ylab = 14, fontsize_axis_text = 11, nbreaks = 2, scientific = FALSE, show_axis, verbose = TRUE){

  if(log10_transform & !is.null(col_mutation_type)){
    if (verbose) cli::cli_alert_warning(
        "{.strong TMB plot}: Ignoring `col_mutation_type` since `log10_transform = TRUE`.
        This is because you cannot accurately plot stacked bars on a logarithmic scale")
    col_mutation_type <- NULL
  }


  if(is.null(col_mutation_type)){
    data[["MutationType"]] <- NA
  }
  else {
    data <- dplyr::rename(data, "MutationType" = {{col_mutation_type}})
  }

  df_counts <- data |>
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
    ggplot2::ggplot(ggplot2::aes(y = Mutations, x = .data[[col_samples]])) +
    ggiraph::geom_col_interactive(
      ggplot2::aes(
        tooltip = .data[["Tooltip"]],
        data_id = .data[[col_samples]],
        fill = .data[["MutationType"]]
      ),
      width = 1,
      show.legend = FALSE
    )

  # Fill palette
  gg <- gg + ggplot2::scale_fill_manual(values = palette, na.value = colour_mutation_type_unspecified)

  #
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
  #trans = "log10"
  labels = ifelse(scientific, yes = scales::label_scientific(), no = scales::label_comma())

  gg <- gg + ggplot2::scale_y_continuous(
    trans = trans,
    oob = scales::oob_squish_any,
    n.breaks = nbreaks,
    labels = labels,
    expand = ggplot2::expansion(c(0, 0))
    )

  # Y Axis Title
  ylabel = ifelse(log10_transform, yes = "log10\nnMuts", no = "nMuts")
  gg <- gg + ggplot2::ylab(ylabel)

  if(!show_ylab)
    gg <- gg + ggplot2::theme(axis.title.y = ggplot2::element_blank())



  # Show/hide axes
  if(!show_axis)
    gg <- gg + ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank()
      )

  return(gg)
}


#' Combine margin plots with main plot
#'
#' @param gg_main main oncoplot tileplot (ggplot)
#' @param gg_tmb barplot describing total mutations. Set to NULL to not draw barplot (ggplot)
#' @param gg_gene barplot describing number of mutated samples per gene. Set to NULL to not draw barplot (ggplot)
#' @param gg_metadata tile plot describing sample-level metadata
#' @param gg_tmb_height percentage of plot height taken up by TMB plot (should be between 5-95) (number)
#' @param gg_gene_width percentage of plot width taken up by genebar plot (should be between 5-95) (number)
#' @param gg_metadata_height percentage of plot height taken up by metadata plot (should be between 5-95) (number)
#' @return patchwork object (or ggplot obj if both `gg_tmb` and `gg_gene` are NULL)
#'
combine_plots <- function(gg_main, gg_tmb = NULL, gg_gene = NULL, gg_metadata = NULL, gg_tmb_height, gg_gene_width, gg_metadata_height){
  assertions::assert(gg_tmb_height + gg_metadata_height < 95)
  assertions::assert(gg_gene_width < 95)

  gg_main_height = 100 - gg_tmb_height - gg_metadata_height
  gg_main_top = gg_tmb_height + 1
  gg_main_bottom = gg_tmb_height + 1 + gg_main_height

  gg_main_width = 100 - gg_gene_width

  #browser()
  # Remove all legends from the cowplot object
  # # Can just comment if statement out plus change metadata plot create return_gglist argument to FALSE to remove cowplot
  # if(is.null(gg_gene) & is.null(gg_tmb) &  !is.null(gg_metadata)){
  #
  #   ls <- gg_metadata$plotlist
  #   ncols_metadata = length(gg_metadata$plotlist)
  #   metadata_height = gg_metadata_height / ncols_metadata
  #   ls[['main']] <- gg_main + ggplot2::theme(legend.position = 'none')
  #   ls <- ls[c('main', names(ls)[names(ls) != 'main'])]
  #   cow = cowplot::plot_grid(plotlist = ls, align = "v", axis = "lr", ncol=1, rel_heights = c(40, rep(5, times=ncols_metadata)))
  #
  #   #browser()
  #   return(cow)
  # }

  # Define layouts (will need to edit to make layout respect gg_main_height, gg_main_width and gg_metadata_height)
  layout <- c(
    patchwork::area(t = gg_main_top, l = 0, b = gg_main_bottom, r = gg_main_width), # Main Plot
    if(!is.null(gg_tmb)) patchwork::area(t = 0, l = 0, b = gg_tmb_height, r = gg_main_width) else patchwork::area(), # TMB Barplot
    if(!is.null(gg_gene)) patchwork::area(t = gg_main_top, l = gg_main_width + 1, b =  gg_main_bottom, r = gg_main_width + gg_gene_width + 1) else patchwork::area(), # Genbar
    if(!is.null(gg_metadata)) patchwork::area(t = gg_main_bottom + 1, l = 0, b = gg_main_bottom + 1 + gg_metadata_height, r = gg_main_width) else patchwork::area() # Metadata
    )

  # Adjust margins of main plot
  gg_main_margins <- gg_main$theme$plot.margin
  unit <- unique(grid::unitType(gg_main_margins))

  gg_main <- gg_main + ggplot2::theme(plot.margin = ggplot2::margin(
    t = ifelse(!is.null(gg_tmb), yes = 0, no = gg_main_margins[1]),
    r = ifelse(!is.null(gg_gene), yes = 0, no = gg_main_margins[2]),
    b = gg_main_margins[3],
    l = gg_main_margins[4],
    unit = unit
  ))

  # Compose final plot
  gg_final <- gg_main + gg_tmb + gg_gene + gg_metadata + patchwork::plot_layout(design = layout, guides = "collect") &
    ggplot2::theme(legend.margin = ggplot2::margin(0, 0, 0 ,0))

  # # Both TMB and gene plots supplied
  # if(!is.null(gg_tmb) & !is.null(gg_gene)){
  #   gg_final <- gg_tmb + patchwork::plot_spacer() + gg_main + gg_gene +
  #     patchwork::plot_layout(
  #       ncol = 2,
  #       widths = c(gg_main_width, gg_gene_width),
  #       heights = c(gg_tmb_height, gg_main_height)
  #     )
  # }
  # # Only TMB
  # else if(!is.null(gg_tmb) & is.null(gg_gene)){
  #   gg_final <- gg_tmb / gg_main +
  #     patchwork::plot_layout(
  #       heights = c(gg_tmb_height, gg_main_height)
  #     )
  # }
  # # Only Gene
  # else if(is.null(gg_tmb) & !is.null(gg_gene)){
  #   gg_final <- gg_main + gg_gene +
  #     patchwork::plot_layout(
  #       ncol = 2,
  #       widths = c(gg_main_width, gg_gene_width)
  #     )
  # }
  # # Neither TMB nor Gene
  # else if(is.null(gg_tmb) & is.null(gg_gene)){
  #   gg_final <- gg_main
  # }
  # else
  #   cli::cli_abort("unexplained case when combining margin plots, package maintainer should please explicitly describe how plots should combine")


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
#' @param data dataframe with a column describing sample IDs (data.frame)
#' @param col_samples name of column in `data` containing sample IDs (character)
#' @param samples_to_show the samples we want to show in plots.
#' These samples should be the only ones represented in data.frame content,
#'  and any missing ones will be added as factor levels (character)
#'
#' @return data.frame
#'
unify_samples <- function(data, col_samples, samples_to_show){
  if(is.null(data)) return(data)

  # Filter to include ONLY samples in samples_to_show
  data <- data[data[[col_samples]] %in% samples_to_show,]

  # Drop any extra levels based on original content
  data[[col_samples]] <- droplevels(as.factor(data[[col_samples]]))

  # add levels for any samples_to_show that are missing from content
  data[[col_samples]] <- forcats::fct_expand(data[[col_samples]], samples_to_show)

  # Ensure metadata columns are in the same order as the sequence of samples_to_show
  data[[col_samples]] <- forcats::fct_relevel(data[[col_samples]], samples_to_show)

  return(data)
}

#' Make strings prettier for printing
#'
#' Takes an input string and 'beautify' by converting underscores to spaces and
#'
#' @param string input string
#'
#' @return string
#'
beautify <- function(string){
  # underscores to spaces
  string <- gsub(x=string, pattern = "_", replacement = " ")

  # camelCase to camel Case
  string <- gsub(x=string, pattern = "([a-z])([A-Z])", replacement = "\\1 \\2")


  # Capitalise Each Word
  string <- gsub(x=string, pattern = "^([a-z])",  perl = TRUE, replacement = ("\\U\\1"))
  string <- gsub(x=string, pattern = " ([a-z])",  perl = TRUE, replacement = (" \\U\\1"))

  return(string)
}

get_genes_for_oncoplot <- function(data, pathway_df = NULL, col_samples, col_genes, topn, genes_to_ignore = NULL, return_extra_genes_if_tied = FALSE, genes_to_include = NULL, verbose = TRUE){
  # Look exclusively at a custom set of genes
  if (!is.null(genes_to_include)) {
    genes_not_found <- genes_to_include[!genes_to_include %in% data[[col_genes]]]

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
      data = data,
      col_samples = col_samples,
      col_genes = col_genes,
      topn = topn,
      return_extra_genes_if_tied = return_extra_genes_if_tied,
      genes_to_ignore = genes_to_ignore,
      verbose = verbose
    )
  }

  # Potentially add pathway sort

  return(genes_for_oncoplot)
}

#' Identify top genes from a mutation df
#'
#' Identify top genes from a mutation df
#'
#' @inheritParams ggoncoplot
#'
#' @return vector of topn genes. Their order will be their rank (most mutated = first) (character)
#'
identify_topn_genes <- function(data, col_samples, col_genes, topn, genes_to_ignore = NULL, return_extra_genes_if_tied = FALSE, verbose = TRUE){
  assertions::assert_flag(return_extra_genes_if_tied)
  if(!is.null(genes_to_ignore)) assertions::assert_character(genes_to_ignore)
  assertions::assert_number(topn)
  assertions::assert_greater_than(topn, minimum = 0)
  assertions::assert_flag(verbose)

  # Identify top genes by frequency
  df_data_gene_counts <- data |>
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
  assertions::assert(is.character(mutated_genes) | is.factor(mutated_genes))
  assertions::assert_character(genes_informing_score)
  assertions::assert_numeric(gene_rank)
  assertions::assert_equal(length(genes_informing_score), length(gene_rank))


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
#' @inheritParams ggoncoplot_options
#' @importFrom ggplot2 %+replace%
theme_oncoplot_default <- function(show_legend_titles = TRUE, fontsize_legend_title = 12,
                                   fontsize_legend_text = 12, ...) {

  ggplot2::theme_bw(...) %+replace%
    ggplot2::theme(
      panel.border = ggplot2::element_rect(linewidth = 1, fill = NA),
      #panel.grid.minor = ggplot2::element_line(colour = "red"),
      panel.grid.major = ggplot2::element_blank(),
      # panel.grid.minor.y = ggplot2::element_line(colour = "red"),
      axis.title = ggplot2::element_text(face = "bold"),
      legend.title = if(show_legend_titles)
        ggplot2::element_text(size = fontsize_legend_title, face = "bold")
      else
        ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = fontsize_legend_text),
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
  assertions::assert_character(colnames)
  assertions::assert_dataframe(data)


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
#' Takes same data input as ggoncoplot and returns a dataframe with 'Sample' and 'Gene' columns
#' ONLY for sample-gene pairs that are unmutated. This lets us colour render them separately (as grey)
#'
#' @inheritParams ggoncoplot_plot
#'
#' @return  a dataframe with 'Sample' and 'Gene' columns ONLY for sample-gene pairs that are unmutated. This lets us colour render them separately (as grey)  (data.frame)
get_nonmutated_tiles <- function(data){
  samples <- levels(data[['Sample']])

  non_mutated_tiles_df  <- expand.grid(
    Sample = samples,
    Gene = unique(data[["Gene"]])
  )

  nomutations <- ! paste0(
    non_mutated_tiles_df[['Sample']],
    non_mutated_tiles_df[['Gene']]
  ) %in%
    unique(
      paste0(
        data[['Sample']],
        data[["Gene"]]
      )
    )

  non_mutated_tiles_df <- non_mutated_tiles_df[nomutations,]

  if(!is.null(data[["Pathway"]])){ # Add pathway col back in so faceting works
    non_mutated_tiles_df[["Pathway"]] <- data[["Pathway"]][match(non_mutated_tiles_df[["Gene"]], data[["Gene"]])]
  }

  return(non_mutated_tiles_df)
}


#' Calculate Pathway-informed Genes Rankings
#'
#' Which genes should appear at the top of the oncoplot?
#' This function takes pathway and gene ranks and returns a list of genes sorted first by pathway then by gene rank.
#' Gene & pathway rankings can be calculated upstream. By default will use their order in gene_pathway_map.
#'
#'
#'
#' @param generanks gene names in the order they should be ranked, where earlier in vector = further up in oncoplot. (character)
#' @param pathwayranks pathway names in the order they should be ranked, where earlier in vector = further up in oncoplot (character)
#' @param gene_pathway_map dataframe where column 1 = gene names and column 2 = pathway names
#' @return gene names, sorted based on order they should appear in oncoplot (first = top). Only returns genes present in generanks (character)
#'
#'
rank_genes_based_on_pathways <- function(gene_pathway_map,
                                         generanks = unique(as.character(gene_pathway_map[[1]])),
                                         pathwayranks = unique(as.character(gene_pathway_map[[2]]))

    ){

    df_pathway <- gene_pathway_map[gene_pathway_map[[1]] %in% generanks,]

    df <- data.frame(Gene = generanks, GeneRanks = seq_along(generanks))
    df[["Pathway"]] <- gene_pathway_map[[2]][match(df$Gene, gene_pathway_map[[1]])]
    df[["PathwayRanks"]] <- match(df$Pathway, pathwayranks)

    df <- df[order(df[['PathwayRanks']], df[['GeneRanks']]), , drop = FALSE]

    return(df[['Gene']])
}

as_pct <- function(x, digits = 1, sep="", multiply_by_100 = TRUE){
  if(multiply_by_100)
    x <- x*100

  x <- round(x, digits = digits)
  paste0(x,sep, "%")
}

# Visual Options ----------------------------------------------------------

#' ggoncoplot options
#'
#' Customise the look of your oncoplot.
#'
#' @param interactive_svg_width dimensions of interactive plot (number)
#' @param interactive_svg_height dimensions of interactive plot (number)
#' @param show_sample_ids show sample_ids_on_x_axis (flag)
#' @param colour_pathway_text colour of text describing pathways (string)
#' @param colour_pathway_bg background fill colour of pathway strips (string)
#' @param colour_pathway_outline outline colour of pathway strips (string)
#' @param pathway_text_angle angle of pathway text label (typically 0 or 90 degrees) (number)
#' @param xlab_title x axis label. Set `xlab_title = NULL` to remove title (string)
#' @param ylab_title y axis of interactive plot. Set `ylab_title = NULL` to remove title (string)
#' @param fontsize_xlab size of x axis title (number)
#' @param fontsize_ylab size of y axis title (number)
#' @param fontsize_genes size of y axis text (gene names) (number)
#' @param fontsize_samples size of x axis text (sample names). Ignored unless show_sample_ids is set to true (number)
#' @param fontsize_count fontsize of gene mutation count x axis (number)
#' @param fontsize_tmb_title fontsize of y axis title for TMB marginal plot (number)
#' @param fontsize_tmb_axis fontsize of y axis text for TMB marginal plot (number)
#' @param fontsize_pathway fontsize of y axis strip text describing gene pathways (number)
#' @param fontsize_metadata_text fontsize of the y axis text for in the sample metadata plot (number)
#' @param fontsize_legend_title fontsize of the legend titles  (number)
#' @param fontsize_legend_text fontsize of the legend text (number)
#' @param fontsize_metadata_legend_text fontsize of the text in metadata legends. Will default to \code{fontsize_legend_title} (number)
#' @param fontsize_metadata_legend_title fontsize of the titles of metadata legends. Will default to \code{fontsize_legend_text} (number)
#' @param fontsize_metadata_barplot_y_numbers fontsize of the text describing numeric barplot max & min values (number)
#' @param tile_height proportion of available vertical space each tile will take up (0-1) (number)
#' @param tile_width proportion of available horizontal space each tile take up (0-1) (number)
#' @param colour_backround colour used for background non-mutated tiles (string)
#' @param colour_mutation_type_unspecified colour of mutations in oncoplot and margin plots if `col_mutation_type` is not supplied (string)
#' @param show_legend show the oncoplot legend
#' @param show_ylab_title show y axis title of oncoplot (flag)
#' @param show_xlab_title show x axis title of oncoplot (flag)
#' @param show_ylab_title_tmb show y axis title of TMB margin plot (flag)
#' @param show_axis_gene show x axis line/ticks/labels for gene barplot (flag)
#' @param show_axis_tmb show y axis line/ticks/labels for TMB barplot (flag)
#' @param log10_transform_tmb log10 transform total number of mutations for TMB marginal plot (flag)
#' @param scientific_tmb display tmb counts in scientific notation (flag)
#' @param plotsize_tmb_rel_height percentage of vertical space TMB margin plot should take up. Must be some value between 5-90 (number)
#' @param plotsize_gene_rel_width percentage of horizontal space the gene barplot should take up. Must be some value between 5-90 (number)
#' @param plotsize_metadata_rel_height percentage of vertical space the metadata tile plot should take up. Must be some value between 5-90 (number)
#' @param ggoncoplot_guide_ncol how many columns to use when describing oncoplot legend (number)
#' @param genebar_label_padding how much padding to add to the x axis of the gene barplot (number)
#' @param genebar_label_nudge how much padding to add between the gene barplot and bar annotations (number)
#' @param genebar_label_round how many digits to round the genebar labels to (number)
#' @param show_genebar_labels should gene barplot be labelled with % of samples the gene is mutated in (flag)
#' @param genebar_only_pad_when_labels_shown only apply \code{genebar_label_padding} when labels are shown (flag)
#' @param genebar_scale_breaks fine-grained control over the x axis breaks on the gene barplot.
#' One of:
#'   - `NULL` for no minor breaks
#'   - `waiver()` for the default breaks (none for discrete, one minor break
#'     between each major break for continuous)
#'   - A numeric vector of positions
#'   - A function that given the limits returns a vector of minor breaks. When
#'     the function has two arguments, it will be given the limits and major
#'     break positions.
#' @param genebar_scale_n_breaks an integer guiding the number of breaks The algorithm
#'   may choose a slightly different number to ensure nice break labels. Will
#'   only have an effect if `genebar_scale_breaks = ggplot2::waiver()`. Use `NULL` to use the default
#' @param show_legend_titles show legend titles (flag)
#' @param metadata_legend_nrow number of rows allowed per metadata legend (number)
#' @param metadata_legend_ncol number of columns allowed per metadata legend (number)
#' @param legend_key_size width of the legend key block (number)
#' @param metadata_legend_key_size width of the legend key block (number). Defaults to \code{legend_key_size}
#' @param metadata_na_marker character used to indicate data is missing (string)
#' @param metadata_na_marker_size size of character used when data is missing (number)
#' @param metadata_maxlevels or categorical variables, what is the maximum number of distinct values to allow (too many will make it hard to find a palette that suits) (number)
#' @param metadata_numeric_plot_type visual representation of numeric properties. One of 'bar', for bar charts, or 'heatmap' for heatmaps
#' @param metadata_legend_orientation_heatmap the orientation of heatmaps in legends. One of "horizontal" or "vertical"
#'   number of breaks given by the transformation.
#' @return ggoncoplot options object ready to be passed to [ggoncoplot()] \code{options} argument
#' @export
#'
#' @examples
#'
#' # Read GBM MAF file
#' gbm_csv <- system.file(
#'   package = "ggoncoplot",
#'   "testdata/GBM_tcgamutations_mc3_maf.csv.gz"
#' )
#' gbm_df <- read.csv(file = gbm_csv, header = TRUE)
#'
#' # Plot Oncoplot and Customise Options
#' gbm_df |>
#'   ggoncoplot(
#'     col_genes = 'Hugo_Symbol',
#'     col_samples = 'Tumor_Sample_Barcode',
#'     col_mutation_type = 'Variant_Classification',
#'
#'     # Customise Visual Options
#'     options = ggoncoplot_options(
#'
#'       # Interactive Plot Options
#'       interactive_svg_width = 12,
#'       interactive_svg_height = 6,
#'
#'       # Relative height of different plotsizes
#'       plotsize_tmb_rel_height = 10,
#'       plotsize_gene_rel_width = 20,
#'       plotsize_metadata_rel_height = 20,
#'
#'       # Axis Titles
#'       xlab_title = "Glioblastoma Samples",
#'       ylab_title = "Top 10 mutated genes",
#'
#'       # Fontsizes
#'       fontsize_xlab = 40,
#'       fontsize_ylab = 40,
#'       fontsize_genes = 16,
#'       fontsize_samples = 12,
#'       fontsize_count = 14,
#'       fontsize_tmb_title = 14,
#'       fontsize_tmb_axis = 11,
#'       fontsize_pathway = 16,
#'
#'       # Customise Tiles
#'       tile_height = 1,
#'       tile_width = 1,
#'       colour_backround = "grey90",
#'       colour_mutation_type_unspecified = "grey10",
#'
#'       # Show different elements
#'       show_sample_ids = FALSE,
#'       show_ylab_title = FALSE,
#'       show_xlab_title = FALSE,
#'       show_ylab_title_tmb = FALSE,
#'       show_axis_gene = TRUE,
#'       show_axis_tmb = TRUE,
#'
#'       # Transformation and label scales
#'       log10_transform_tmb = TRUE,
#'       scientific_tmb = FALSE,
#'
#'       # Gene Barplot Specific Options
#'       show_genebar_labels = TRUE,
#'       genebar_label_padding = 0.2,
#'       genebar_only_pad_when_labels_shown = TRUE,
#'       genebar_label_nudge = 2,
#'       genebar_label_round = 1,
#'
#'       # Pathway Faceting Colours / Text
#'       colour_pathway_text = "white",
#'       colour_pathway_bg = "grey10",
#'       colour_pathway_outline = "black",
#'       pathway_text_angle = 0,
#'
#'       # Legend number of columns
#'       ggoncoplot_guide_ncol = 2
#'     )
#'   )
ggoncoplot_options <- function(
    # Interactive Plot Options
    interactive_svg_width = 12,
    interactive_svg_height = 6,

    # Relative height of different plotsizes
    plotsize_tmb_rel_height = 10,
    plotsize_gene_rel_width = 20,
    plotsize_metadata_rel_height = 20,

    # Axis Titles
    xlab_title = "Sample",
    ylab_title = "Gene",

    # Fontsizes
    fontsize_xlab = 26,
    fontsize_ylab = 26,
    fontsize_genes = 16,
    fontsize_samples = 12,
    fontsize_count = 14,
    fontsize_tmb_title = 14,
    fontsize_tmb_axis = 11,
    fontsize_pathway = 16,
    fontsize_legend_title = 12,
    fontsize_legend_text = 12,

    # Customise Tiles
    tile_height = 1,
    tile_width = 1,
    colour_backround = "grey90",
    colour_mutation_type_unspecified = "grey10",

    # Show different elements
    show_sample_ids = FALSE,
    show_ylab_title = FALSE,
    show_xlab_title = FALSE,
    show_ylab_title_tmb = FALSE,
    show_legend = TRUE,
    show_legend_titles = TRUE,
    show_axis_gene = TRUE,
    show_genebar_labels = FALSE,
    show_axis_tmb = TRUE,

    # Transformation and label scales
    log10_transform_tmb = TRUE,
    scientific_tmb = FALSE,

    # Gene Barplot Specific Options
    genebar_label_padding = 0.3,
    genebar_only_pad_when_labels_shown = TRUE,
    genebar_label_nudge = 2,
    genebar_label_round = 0,
    genebar_scale_breaks = ggplot2::waiver(),
    genebar_scale_n_breaks = 3,

    # Pathway Faceting Colours / Text
    colour_pathway_text = "white",
    colour_pathway_bg = "grey10",
    colour_pathway_outline = "black",
    pathway_text_angle = 0,

    # Legend number of columns
    ggoncoplot_guide_ncol = 2,
    legend_key_size = 0.4,

    # ====== Metadata ======
    # Metadata: Fontsizes
    fontsize_metadata_text = 12, # Y axis text
    fontsize_metadata_legend_title = fontsize_legend_title,
    fontsize_metadata_legend_text = fontsize_legend_text,
    fontsize_metadata_barplot_y_numbers = 8,

    # Metadata: Legend Layout
    metadata_legend_nrow = NULL,
    metadata_legend_ncol = NULL,
    metadata_legend_key_size = legend_key_size,

    # Metadata: Dealing with NAs
    metadata_na_marker = "!",
    metadata_na_marker_size = 8,

    # Metadata: Processing Levels
    metadata_maxlevels = 6,

    # Metadata: Numeric Values
    metadata_numeric_plot_type = c("bar", "heatmap"),
    metadata_legend_orientation_heatmap = c("horizontal", "vertical")
  ){


  # Assertions --------------------------------------------------------------
  assertions::assert_number(fontsize_xlab)
  assertions::assert_number(fontsize_ylab)
  assertions::assert_number(fontsize_genes)
  assertions::assert_number(fontsize_samples)
  assertions::assert_number(fontsize_tmb_title)
  assertions::assert_number(fontsize_tmb_axis)
  assertions::assert_number(tile_height)
  assertions::assert_number(tile_width)
  assertions::assert_string(colour_backround)
  assertions::assert_number(fontsize_legend_title)
  assertions::assert_number(fontsize_legend_text)
  assertions::assert_flag(log10_transform_tmb)
  assertions::assert_string(colour_mutation_type_unspecified)
  assertions::assert_flag(scientific_tmb)
  assertions::assert(dplyr::between(plotsize_gene_rel_width, 5, 90), msg = "plotsize_gene_rel_width must be between 5 & 90 (inclusive).")
  assertions::assert(dplyr::between(plotsize_tmb_rel_height, 5, 90), msg = "plotsize_tmb_rel_height must be between 5 & 90 (inclusive).")
  assertions::assert_flag(show_axis_gene)
  assertions::assert_flag(show_axis_tmb)
  assertions::assert_flag(show_genebar_labels)
  assertions::assert_number(genebar_label_padding)
  assertions::assert_number(genebar_label_nudge)
  assertions::assert_flag(genebar_only_pad_when_labels_shown)
  assertions::assert_number(genebar_label_round)
  assertions::assert_greater_than_or_equal_to(genebar_label_round, minimum = 0)
  assertions::assert_flag(show_legend)
  assertions::assert_flag(show_legend_titles)
  assertions::assert_number(legend_key_size)
  assertions::assert_number(fontsize_metadata_text)

  # Metadata options
  if(!is.null(fontsize_metadata_legend_title)) assertions::assert_number(fontsize_metadata_legend_title)
  if(!is.null(fontsize_metadata_legend_text)) assertions::assert_number(fontsize_metadata_legend_text)
  if(!is.null(metadata_legend_nrow)) assertions::assert_whole_number(metadata_legend_nrow)
  if(!is.null(metadata_legend_ncol)) assertions::assert_whole_number(metadata_legend_ncol)
  assertions::assert_numeric(fontsize_metadata_barplot_y_numbers)
  assertions::assert_number(metadata_legend_key_size)
  assertions::assert_number(metadata_na_marker_size)
  assertions::assert_whole_number(metadata_maxlevels)
  metadata_numeric_plot_type <- rlang::arg_match(metadata_numeric_plot_type)
  metadata_legend_orientation_heatmap <- rlang::arg_match(metadata_legend_orientation_heatmap)


  options <- list(
    interactive_svg_width = interactive_svg_width,
    interactive_svg_height = interactive_svg_height,
    plotsize_tmb_rel_height = plotsize_tmb_rel_height,
    plotsize_gene_rel_width = plotsize_gene_rel_width,
    plotsize_metadata_rel_height = plotsize_metadata_rel_height,
    xlab_title = xlab_title,
    ylab_title = ylab_title,
    fontsize_xlab = fontsize_xlab,
    fontsize_ylab = fontsize_ylab,
    fontsize_genes = fontsize_genes,
    fontsize_samples = fontsize_samples,
    fontsize_count = fontsize_count,
    fontsize_tmb_title = fontsize_tmb_title,
    fontsize_tmb_axis = fontsize_tmb_axis,
    fontsize_pathway = fontsize_pathway,
    fontsize_legend_title = fontsize_legend_title,
    fontsize_legend_text = fontsize_legend_text,
    fontsize_metadata_text = fontsize_metadata_text,
    tile_height = tile_height,
    tile_width = tile_width,
    colour_backround = colour_backround,
    colour_mutation_type_unspecified = colour_mutation_type_unspecified,
    show_sample_ids = show_sample_ids,
    show_legend = show_legend,
    show_ylab_title = show_ylab_title,
    show_xlab_title = show_xlab_title,
    show_ylab_title_tmb = show_ylab_title_tmb,
    show_axis_gene = show_axis_gene,
    show_genebar_labels = show_genebar_labels,
    show_axis_tmb = show_axis_tmb,
    log10_transform_tmb = log10_transform_tmb,
    scientific_tmb = scientific_tmb,
    colour_pathway_text = colour_pathway_text,
    colour_pathway_bg = colour_pathway_bg,
    colour_pathway_outline = colour_pathway_outline,
    pathway_text_angle = pathway_text_angle,
    ggoncoplot_guide_ncol = ggoncoplot_guide_ncol,
    genebar_label_padding = genebar_label_padding,
    genebar_only_pad_when_labels_shown = genebar_only_pad_when_labels_shown,
    genebar_label_nudge = genebar_label_nudge,
    genebar_label_round = genebar_label_round,
    genebar_scale_breaks = genebar_scale_breaks,
    genebar_scale_n_breaks = genebar_scale_n_breaks,
    show_legend_titles = show_legend_titles,
    legend_key_size = legend_key_size,
    fontsize_metadata_legend_title = fontsize_metadata_legend_title,
    fontsize_metadata_legend_text = fontsize_metadata_legend_text,
    metadata_legend_nrow = metadata_legend_nrow,
    metadata_legend_ncol = metadata_legend_ncol,
    fontsize_metadata_barplot_y_numbers = fontsize_metadata_barplot_y_numbers,
    metadata_legend_key_size = metadata_legend_key_size,
    metadata_na_marker_size = metadata_na_marker_size,
    metadata_maxlevels = metadata_maxlevels,
    metadata_numeric_plot_type = metadata_numeric_plot_type,
    metadata_legend_orientation_heatmap = metadata_legend_orientation_heatmap
  )

  class(options) <- "ggoncoplot_options"

  return(options)
}
