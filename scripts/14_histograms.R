# Histograms of cluster size vs. number of clusters
# Brendan Furneaux
# May 2021

hist_plan <- list(
  #### hist_plot ####
  hist_plot =tar_target(
    hist_plot,
    cluster_key %>%
      # dplyr::select(-seq_id, -reads_type) %>%
      dplyr::mutate(cluster_type = sub("_.*", "", OTU)) %>%
      tidyr::pivot_longer(cols = GH90:SH99, names_to = "sh_type", values_to = "SH") %>%
      dplyr::group_by(cluster_type, sh_type, SH) %>%
      dplyr::summarize(OTUs = dplyr::n_distinct(OTU)) %>%
      recode_cluster_types() %>%
      dplyr::mutate_at("sh_type", factor,
                       levels = c("SH99", "SH97", "GH90")) %>%
      ggplot(aes(x = OTUs, fill = cluster_type)) +
      facet_grid(sh_type ~ cluster_type, scales = "free_y") +
      geom_bar(width = 0.03) +
      scale_y_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000),
                         minor_breaks = sqrt(10) * c(1, 10, 100, 1000),
                         name = "Number of clusters") +
      scale_x_log10(name = "Number of OTUs per cluster",
                    breaks = c(1, 3, 10, 30, 100)) +
      scale_fill_discrete(guide = NULL) +
      theme_bw() +
      theme(strip.background = element_blank()),
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
        device = fun, width = 4, height = 3.5, dpi = 150
      )
    )
  )
)
