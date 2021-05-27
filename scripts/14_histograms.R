# Histograms of cluster size vs. number of clusters
# Brendan Furneaux
# May 2021

hist_plan <- list(
  #### hist_plot ####
  hist_plot =tar_target(
    hist_plot,
    dplyr::left_join(
      reads_tab %>%
        tidyr::pivot_longer(cols = c("ampliseq", "vsearch", "single_link"),
                            names_to = "cluster_type", values_to = "reads"),
      dplyr::select(cluster_key, -OTU) %>%
        unique(),
      by = c("seq_id" = "ITS2_hash")) %>%
      # dplyr::select(-seq_id, -reads_type) %>%
      tidyr::pivot_longer(cols = GH90:SH99, names_to = "sh_type", values_to = "SH") %>%
      dplyr::group_by(cluster_type, sh_type, SH) %>%
      dplyr::summarize(reads = sum(reads)) %>%
      recode_cluster_types() %>%
      dplyr::mutate_at("cluster_type", forcats::fct_relabel, chartr,
                       old = " ", new = "\n") %>%
      dplyr::mutate_at("sh_type", factor,
                       levels = c("GH90", "SH97", "SH99", "OTU")) %>%
      ggplot(aes(x = log10(reads), weight = reads)) +
      facet_grid(cluster_type ~ sh_type, scales = "free_y") +
      geom_histogram(bins = 20) +
      scale_x_continuous(breaks = 0:4,
                         labels = c("1", "10", "100", "1k", "10k"),
                         name = "reads per cluster") +
      scale_y_continuous(name = "total reads", breaks = 2000 * 0:4,
                         labels = c("0", paste0(2 * 1:4, "k"))),
    packages = "ggplot2"
  ),
  tar_map(
    values = plot_type_meta,
    names = ext,
    #### hist_plot_file_{ext} ####
    hist_plot_file = tar_file(
      hist_plot_file,
      write_and_return_file(
        hist_plot,
        file.path(figdir, sprintf("histograms.%s", ext)),
        device = fun, width = 6.25, height = 4, dpi = 150
      )
    )
  )
)
