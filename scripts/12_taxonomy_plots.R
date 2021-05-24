library(targets)
library(tarchetypes)
library(magrittr)

phylotax_meta <-
  tibble::tibble(
    group = c("protists", "fungi"),
    rank = c("kingdom", "phylum"),
    tree = c("tree_protists", "tree_fungi_new"),
    cutoff = c(0.015, 0.015)
  ) %>%
  dplyr::mutate_at("tree", rlang::syms)

taxplot_meta <-
  dplyr::left_join(
    dplyr::mutate(
      phylotax_meta,
      tax_consensus = rlang::syms(paste0("tax_consensus_", group))
    ),
    tibble::tibble(
      clust_type = c("ampliseq", "swarm", "vsearch"),
      cl_id = c("otuA", "otuS", "otuC"),
      regions = c("regions_as", "regions_sl", "regions_vs")
    ),
    by = character()
  ) %>%
  dplyr::mutate(
    otutab = sprintf("otu_table_%s_%s", clust_type, group)
  ) %>%
  dplyr::mutate_at(c("tree", "otutab", "rank", "regions"), rlang::syms)

phylotax_plan <- tar_map(
  values = phylotax_meta,
  names = group,
  # make phylogenetic consensus assignments
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
  )
)

taxplot_plan <- tar_map(
  values = taxplot_meta,
  names = c(group, cl_id),
  # data for taxonomy plot
  tar_fst_tbl(
    taxplot_data,
    # convert OTU table to relative abundances
    otutab %>%
      tibble::column_to_rownames("OTU") %>%
      vegan::decostand(method = "total", MARGIN = 2) %>%
      tibble::as_tibble(rownames = "OTU") %>%
      # convert to long format
      # i.e. each OTU/sample combination on its own row (removing empties)
      tidyr::pivot_longer(-1, names_to = "sample", values_to = "reads") %>%
      dplyr::filter(reads > 0) %>%
      # add tree tip label (5.8S_LSU) for each OTU
      dplyr::left_join(
        dplyr::select(regions, OTU = "seq_id", "label"),
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
        file.path(figdir, sprintf("taxonomy_%s_%s.%s", cl_id, group, ext)),
        device = fun, width = 6.25, height = 4, dpi = 150
      )
    )
  ),
  tar_qs(
    samples_physeq,
    phyloseq::phyloseq(
      phyloseq::otu_table(
        tibble::column_to_rownames(otutab, "OTU"),
        taxa_are_rows = TRUE
      ),
      phyloseq::tax_table(
        taxplot_data %>%
          dplyr::select(OTU, kingdom:genus) %>%
          unique() %>%
          tibble::column_to_rownames("OTU") %>%
          as.matrix()
      ),
      phyloseq::sample_data(samples_df)
    )
  ),
  tar_file(
    physeq_file,
    write_and_return_file(
      samples_physeq,
      sprintf("output/data/phyloseq_%s_%s.rds", cl_id, group)
    )
  )
)
