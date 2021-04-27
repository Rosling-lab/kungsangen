venn_targets <- list(
  venn_data = tar_fst_tbl(
    venn_data,
    list(
      "fungi" = otu_table_ampliseq_fungi,
      "protists" = otu_table_ampliseq_protists
    ) %>%
      dplyr::bind_rows(.id = "group") %>%
      # convert to long format
      # i.e. each OTU/sample combination on its own row (removing empties)
      tidyr::pivot_longer(-(1:2), names_to = "sample", values_to = "reads") %>%
      # convert to relative abundance
      dplyr::group_by(sample) %>%
      dplyr::mutate(reads = reads/sum(reads)) %>%
      # add sample info
      dplyr::left_join(
        tibble::rownames_to_column(samples_df, "sample"),
        by = "sample"
      ) %>%
      dplyr::group_by(group, OTU, Sites) %>%
      dplyr::summarize(reads = mean(reads)) %>%
      dplyr::mutate(
        category = dplyr::case_when(
          sum(reads > 0 & Sites == "Dry") == 0 ~ "Wet-only",
          sum(reads > 0 & Sites == "Wet") == 0 ~ "Dry-only",
          TRUE ~ "Common"
        ) %>%
          factor(levels = c("Dry-only", "Common", "Wet-only"))
      ) %>%
      dplyr::group_by(group, Sites, category) %>%
      dplyr::summarize(OTUs = sum(reads > 0), reads = sum(reads)) %>%
      dplyr::filter(reads > 0)
  ),
  venn_barplot = tar_target(
    venn_barplot,
    venn_data %>%
      dplyr::mutate(
        Sites = paste(Sites, "sites")#,
        # category = forcats::fct_relabel(category, paste, "OTUs")
      ) %>%
      ggplot(aes(x = category, y = reads, fill = group, label = formatC(reads, digits = 2, format = "f"))) +
      geom_col(position = "stack", width = 0.5) +
      geom_text(position = position_stack(vjust = 0.5)) +
      scale_fill_discrete(name = NULL) +
      facet_wrap(~Sites, scales = "free_x", strip.position = "bottom", ) +
      ylab("Relative abundance") +
      xlab(NULL) +
      theme(strip.placement = "outside",
            strip.background = element_rect(fill = NA, color = "black"),
            legend.position = "bottom"),
    packages = "ggplot2"
  ),
  tar_map(
    values = list(
      g = c("fungi", "protists"),
      color = c("tomato", "cyan")
    ),
    names = g,
    tar_target(
      venn_plot,
      venn_data %>%
        dplyr::filter(group == g) %>%
        dplyr::group_by(category) %>%
        dplyr::summarize(OTUs = max(OTUs)) %$%
        eulerr::euler(
          combinations = c(Dry = OTUs[1], "Dry&Wet" = OTUs[2], Wet = OTUs[3])
        ) %>%
        plot(
          quantities = TRUE,
          fills = list(fill = paste(color, c(1, 4)), alpha = c(0.5, 0.5))
        )
    )
  ),
  tar_target(
    venn_fullplot,
    ggpubr::ggarrange(
      venn_plot_fungi,
      venn_plot_protists,
      venn_barplot,
      ncol = 1,
      heights = c(1, 1, 1.5),
      labels = "AUTO"
    )
  ),
  tar_map(
    values = plot_type_meta,
    names = ext,
    tar_file(
      vennplotfile,
      write_and_return_file(
        venn_fullplot,
        file.path(figdir, sprintf("venn.%s", ext)),
        device = fun, width = 4, height = 7, dpi = 150
      )
    )
  ),
  tar_file(
    venn_xlsx,
    write_and_return_file(venn_data, "output/data/venn.xlsx", "xlsx")
  )
)
