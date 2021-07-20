library(targets)
library(tarchetypes)
library(magrittr)

detection_meta <- tibble::tibble(
  threshold = c(90, 97, 99),
  clusters = paste0("its2_cluster_", threshold)
) %>%
  dplyr::mutate_at("clusters", rlang::syms)

detection_plan <- tar_map(
  values = detection_meta,
  names = threshold,
  #### detection_data_{threshold} ####
  tar_target(
    detection_data,
    parse_clusters(clusters, its2_precluster_singletons) %>%
      dplyr::left_join(reads_tab, by = "seq_id") %>%
      dplyr::group_by(cluster) %>%
      dplyr::summarize(
        n_asvs = sum(ampliseq > 0),
        dplyr::across(where(is.numeric), sum)
      ) %>%
      dplyr::mutate(
        n = pmax(ampliseq, vsearch, single_link),
        threshold = as.character(threshold)
      )
  )
)

detection_plan <- c(
  detection_plan,
  list(
    #### detection_plot ####
    detection_plot = tar_combine(
      detection_plot,
      detection_plan$detection_data,
      command = dplyr::bind_rows(!!!.x) %>%
        ggplot(aes(x = n, y = n_asvs, color = threshold, group = threshold)) +
        geom_point(alpha = 0.2) +
        scale_x_log10() +
        coord_trans(y = "sqrt") +
        scale_color_brewer(type = "qual", palette = 2) +
        stat_smooth(method = glm, method.args = list(family = "poisson")) +
        ylab("ASVs / SH") +
        xlab("Reads / SH"),
      packages = "ggplot2"
    ),
    detection_plotfile = tar_map(
      values = plot_type_meta,
      names = ext,
      #### detection_plotfile_{ext} ####
      tar_file(
        detection_plotfile,
        write_and_return_file(
          detection_plot,
          file.path(figdir, sprintf("detection.%s", ext)),
          device = fun,
          width = 6.25, height = 4, dpi = 150
        )
      )
    )
  )
)
