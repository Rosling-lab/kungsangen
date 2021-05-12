library(targets)
library(tarchetypes)
library(magrittr)

taxplot_meta <- tibble::tibble(
  group = c("protists", "fungi"),
  rank = c("kingdom", "phylum"),
  tree = c("tree_protists", "tree_fungi_new"),
  cutoff = c(0.015, 0.015),
  otu_table = paste0("otu_table_ampliseq_", group)
) %>%
  dplyr::mutate_at(c("tree", "otu_table", "rank"), rlang::syms)

taxplot_plan <- tar_map(
  values = taxplot_meta,
  names = group,
  # make phylogenetic consensus assignments for the fungi
  tar_target(
    tax_consensus,
    phylotax::phylotax(
      tree = tree,
      taxa = tax_all %>%
        dplyr::left_join(
          dplyr::select(hash_key, "full_hash", "label"),
          by = "full_hash"
        ) %>%
        dplyr::filter(label %in% tree$tip.label)
    )
  ),
  # data for taxonomy plot
  tar_fst_tbl(
    taxplot_data,
    # convert OTU table to relative abundances
    otu_table %>%
      tibble::column_to_rownames("OTU") %>%
      vegan::decostand(method = "total", MARGIN = 2) %>%
      tibble::as_tibble(rownames = "OTU") %>%
      # convert to long format
      # i.e. each OTU/sample combination on its own row (removing empties)
      tidyr::pivot_longer(-1, names_to = "sample", values_to = "reads") %>%
      dplyr::filter(reads > 0) %>%
      # add tree tip label (5.8S_LSU) for each OTU
      dplyr::left_join(
        dplyr::select(regions_as, OTU = "seq_id", "label"),
        by = "OTU"
      ) %>%
      # add taxonomy in wide format
      dplyr::inner_join(
        dplyr::select(tax_consensus$assigned, "label", "rank", "taxon") %>%
          tidyr::pivot_wider(names_from = "rank", values_from = "taxon"),
        by = "label"
      ) %>%
      # add sample info
      dplyr::left_join(
        tibble::rownames_to_column(samples_df, "sample"),
        by = "sample"
      )
  ),
  # make plots for reads and OTUs
  tar_map(
    values = list(plottype = rlang::syms(c("reads", "OTUs"))),
    tar_target(
      taxplot,
      taxon_plot(taxplot_data, rank = rank, y = plottype, x = Sites,
                 cutoff = cutoff, cutoff_type = "either"),
      packages = "ggplot2"
    )
  ),
  # combine reads plot and OTUs plot
  tar_target(
    taxplot,
    ggpubr::ggarrange(
      taxplot_reads,
      taxplot_OTUs,
      ncol = 1,
      labels = "AUTO",
      legend = "bottom",
      common.legend = TRUE
    ),
    packages = "ggplot2"
  ),
  tar_map(
    values = plot_type_meta,
    names = ext,
    tar_file(
      taxplotfile,
      write_and_return_file(
        taxplot,
        file.path(figdir, sprintf("taxonomy_%s.%s", group, ext)),
        device = fun, width = 6.25, height = 4, dpi = 150
      )
    )
  )
)
