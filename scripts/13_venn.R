venn_meta <- tibble::tibble(
  group_var_name = c("moisture", "Fritillaria"),
  group_var = rlang::syms(c("Sites", "Sample_type")),
  level1 = c("Dry", "Fritillaria"),
  level1label = c('Dry', 'Fritillaria'),
  level1parselabel = c('"Dry"', 'italic("Fritillaria")'),
  level2 = c("Wet", "Non_Fritillaria"),
  level2label = c('Wet', 'Non-Fritillaria'),
  level2parselabel = c('"Wet"', '"Non-"*italic("Fritillaria")'),
  both = sprintf("%s&%s", level1, level2),
  only1 = paste0(sub("illaria", "", level1label), "-only"),
  only2 = paste0(sub("illaria", "", level2label), "-only")
)

venn_plan <- tar_map(
  values = venn_meta,
  names = group_var_name,
  #### venn_data_{group_var_name} ####
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
      dplyr::group_by(group, OTU, group_var) %>%
      dplyr::summarize(reads = mean(reads)) %>%
      dplyr::mutate(
        category = dplyr::case_when(
          sum(reads > 0 & group_var == level1) == 0 ~ only2,
          sum(reads > 0 & group_var == level2) == 0 ~ only1,
          TRUE ~ "Common"
        ) %>%
          factor(levels = c(only1, "Common", only2))
      ) %>%
      dplyr::group_by(group, group_var, category) %>%
      dplyr::summarize(OTUs = sum(reads > 0), reads = sum(reads)) %>%
      dplyr::filter(reads > 0)
  ),
  #### venn_barplot_{group_var_name} ####
  venn_barplot = tar_target(
    venn_barplot,
    venn_data %>%
      dplyr::mutate(
        facetvar = factor(group_var, levels = c(level1, level2),
                          labels = paste0(c(level1parselabel, level2parselabel), '~"sites"')),
        group = stringr::str_to_title(group),
        # category = forcats::fct_relabel(category, paste, "OTUs")
      ) %>%
      ggplot(aes(x = category, y = reads, fill = group,
                 label = formatC(reads, digits = 2, format = "f"))) +
      geom_col(position = "stack", width = 0.5, color = "white") +
      geom_text(position = position_stack(vjust = 0.5)) +
      scale_fill_manual(values = list(Fungi = "tomato3", Protists = "cyan3"),
                        name = NULL) +
      facet_wrap(~facetvar, scales = "free_x", strip.position = "bottom",
                 labeller = label_parsed) +
      ylab("Relative abundance") +
      xlab(NULL) +
      theme(strip.placement = "outside",
            strip.background = element_rect(fill = NA, color = "black"),
            legend.position = "bottom"),
    packages = "ggplot2",
    tidy_eval = FALSE
  ),
  tar_map(
    values = list(
      g = c("fungi", "protists"),
      color = c("tomato", "cyan")
    ),
    names = g,
    #### venn_plot_{g}_{group_var_name} ####
    tar_target(
      venn_plot,
      venn_data %>%
        dplyr::filter(group == g) %>%
        dplyr::group_by(category) %>%
        dplyr::summarize(OTUs = max(OTUs)) %$%
        eulerr::euler(
          combinations = set_names(OTUs, c(level1, both, level2))
        ) %>%
        plot(
          quantities = TRUE,
          fills = list(fill = paste(color, c(1, 4))),
          edges = list(col = "white", lex = 2),
          labels = c(level1label, level2label)
        ),
      tidy_eval = TRUE
    )
  ),
  #### venn_fullplot_{group_var_name} ####
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
    #### vennplotfile_{ext}_{group_var_name} ####
    tar_file(
      vennplotfile,
      write_and_return_file(
        venn_fullplot,
        file.path(figdir, sprintf("venn_%s.%s", group_var_name, ext)),
        device = fun, width = 4, height = 7, dpi = 150
      )
    )
  ),
  #### venn_xlsx_{group_var_name} ####
  tar_file(
    venn_xlsx,
    write_and_return_file(venn_data, "output/data/venn.xlsx", "xlsx")
  )
)
