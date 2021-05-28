# Generate species accumulation curves
library(magrittr)
library(targets)
library(tarchetypes)

accum_meta <- tibble::tibble(
  clustype = c("as", "vs", "sl"),
  otutab = rlang::syms(c("ampliseq_table", "table_vs", "table_sl")),
  regions = paste0("regions_", clustype)
)

list(
  samples_file = tar_file(samples_file, file.path("start_files", "meta.txt")),
  #### samples_df ####
  # read the metadata
  samples_df = tar_target(
    samples_df,
    read.table(
      samples_file,
      header = TRUE,
      row.names = 1,
      stringsAsFactors = FALSE
    ) %>%
      tidyr::unite("Condition", Sites, Sample_type, sep = "_", remove = FALSE)
  ),

  ### *Alpha diversity at each sample* ####
  #### *begin accum_map {clustype}* ####
  accum_map = tar_map(
    accum_meta,
    #### physeq_sample_{clustype} ####
    # create a phyloseq object
    physeq_sample = tar_target(
      physeq_sample,
      phyloseq::phyloseq(
        phyloseq::otu_table(
          tibble::column_to_rownames(otutab, "OTU"),
          taxa_are_rows = TRUE
        ),
        phyloseq::sample_data(samples_df),
        phyloseq::tax_table(
          dplyr::left_join(dplyr::select(otutab, "OTU"), cluster_key, by = "OTU") %>%
            dplyr::select(GH90, SH97, SH99, ITS2_hash, OTU) %>%
            tibble::column_to_rownames("OTU") %>%
            magrittr::inset2("OTU", value = row.names(.)) %>%
            as.matrix()
        )
      )
    ),
    #### *begin clustmap {accumclust}* ####
    clustmap = tar_map(
      values = list(accumclust = c("OTU", "GH90", "SH97", "SH99")),
      #### clust_physeq_{accumclust}_{clustype} ####
      # make an accumulation curve for reads at each sample
      clust_physeq = tar_target(
        clust_physeq,
        phyloseq::tax_glom(physeq_sample, taxrank = accumclust)
      ),
      #### seq_accum_sample_{accumclust}_{clustype} ####
      seq_accum_sample = tar_target(
        seq_accum_sample,
        # extrapolate the curves up to 5000 reads, with 100 points on each curve
        phyloseq::otu_table(clust_physeq) %>%
          as("matrix") %>%
          iNEXT::iNEXT(endpoint = 5000, knots = 100)
      ),

      #### seq_accum_sample_means_{accumclust}_{clustype} ####
      # take the mean of all the samples in each type
      seq_accum_sample_means = tar_target(
        seq_accum_sample_means,
        accumulation_means(
          seq_accum_sample,
          newx = c(1, 1:100 * 50),
          sampledata = samples_df,
          groupname = "Sites"
        ),
        packages = "iNEXT"
      ),
      #### sample_accum_plot_{accumclust}_{clustype} ####
      sample_accum_plot = tar_target(
        sample_accum_plot,
        # make the plot
        # get the plotting data
        fortify(seq_accum_sample) %>%
          # add forest type annotations
          dplyr::mutate(Type = samples_df[site, "Sites"]) %>%
          # plot it
          ggplot(aes(x = x, y = y)) +
          # lines for each sample
          geom_line(aes(
            group = paste(site, method),
            # darker when interpolated, lighter when extrapolated
            alpha = dplyr::case_when(method == "interpolated" ~ 0.35, TRUE ~0.10)
          )) +
          # use the alpha values I provided as the literal alpha values
          scale_alpha_identity() +
          # small points to represent the actual observed read depth and richness
          geom_point(data = ~dplyr::filter(., method == "observed"), alpha = 0.4, size = 0.5) +
          # thick dashed line for the mean of all the accumulation curves
          geom_line(data = seq_accum_sample_means, size = 0.8, linetype = "dashed") +
          # seperate facet for each forest type.
          # They are side-by-side so we can compare them visually.
          facet_wrap(facets = ~Type, nrow = 1) +
          xlab("reads") +
          scale_x_continuous() +
          ylab("OTU richness"),
        packages = c("ggplot2", "iNEXT")
      ),
      #### *Alpha diversity at each condition* ####

      #### seq_accum_site_{accumclust}_{clustype} ####
      # Calculate for 2 conditions: Wet, Dry
      # make an accumulation curve for reads at each site
      # i.e., what would happen if we had the same samples,
      # but sequenced them more (or less)?
      seq_accum_site = tar_target(
        seq_accum_site,
        phyloseq::merge_samples(clust_physeq, "Sites", fun = dplyr::first) %>%
          phyloseq::otu_table() %>%
          as("matrix") %>%
          t() %>%
          iNEXT::iNEXT(endpoint = 30000, knots = 100)
      ),
      #### samp_accum_site_{accumclust}_{clustype} ####
      # make an accumulation curve for samples at each site
      # i.e., what would happen if we had the same sequencing depth per sample,
      # but got more (or less) samples per site?
      samp_accum_site = tar_target(
        samp_accum_site,
        phyloseq::otu_table(clust_physeq) %>%
          as("matrix") %>%
          vegan::decostand(method = "pa") %>%
          t() %>%
          split.data.frame(samples_df$Sites) %>%
          lapply(t) %>%
          iNEXT::iNEXT(datatype = "incidence_raw", endpoint = 20)
      ),

      #### accum_site_{accumclust}_{clustype} ####
      # make plotting data for both types of accumulation curve at each site
      accum_site = tar_target(
        accum_site,
        list(seq_accum_site, samp_accum_site) %>%
          # get plotting data for each
          lapply(ggplot2::fortify) %>%
          # add a column "xvar" to identify them, and put them together in one tibble
          list(., xvar = c("reads", "samples")) %>%
          purrr::pmap_dfr(tibble::add_column) %>%
          # duplicate the "observed" points as "interpolated" and "extrapolated"
          # so that there isn't a gap in the lines in the plot.
          # this noticable before because they were very closely spaced.
          dplyr::bind_rows(
            dplyr::filter(., method == "observed") %>%
              dplyr::select(-method) %>%
              tidyr::crossing(method = c("interpolated", "extrapolated"))
          ) %>%
          # we will label the observed points. but nothing else
          dplyr::mutate(label = ifelse(method == "observed", as.character(site), ""))
      ),

      #### asymp_site_{accumclust}_{clustype} ####
      # Find the asymptotic species richness estimates
      asymp_site = tar_target(
        asymp_site,
        list(seq_accum_site, samp_accum_site) %>%
          purrr::map(
            ~ cbind(
              dplyr::filter(.$AsyEst, Diversity == "Species richness"),
              .$DataInfo
            ) %>%
              dplyr::select(Estimator, site)
          ) %>%
          # add a column "xvar" to identify them, and put them together in one tibble
          list(., xvar = c("reads", "samples")) %>%
          purrr::pmap_dfr(tibble::add_column)
      ),

      #### *write graphics files* ####
      #### *begin map {ext}* ####
      tar_map(
        values = plot_type_meta,
        #### accumplot1_{ext}_{accumclust}_{clustype} ####
        accumplot1 = tar_file(
          accumplot1,
          file.path(figdir, sprintf("accum1_%s_%s.%s", clustype, accumclust, ext)) %T>%
            ggplot2::ggsave(filename = ., plot = sample_accum_plot, device = fun,
                            width = 6.25, height = 3, dpi = 150)
        ),
        names = ext
      ) #### *end map {ext}* ####
    ), #### *end clustmap {accumclust}* ####

    names = clustype
  ) #### *end accum_map {clustype}* ####
) ->
  accumulation_plan

accumulation_plan <- c(
  accumulation_plan,
  tar_plan(
    #### seq_accum_sample ####
    tar_combine(
      seq_accum_sample,
      accumulation_plan$accum_map[
        startsWith(names(accumulation_plan$accum_map), "seq_accum_sample") &
          !startsWith(names(accumulation_plan$accum_map), "seq_accum_sample_means") &
          !endsWith(names(accumulation_plan$accum_map), "OTU")
      ],
      command = list(!!!.x) %>%
        lapply(ggplot2::fortify) %>%
        dplyr::bind_rows(.id = "cluster_type") %>%
        dplyr::mutate_at("cluster_type", stringr::str_remove, "seq_accum_sample_") %>%
        split_sh_cluster_types(),
      packages = "iNEXT"
    ),
    #### seq_accum_sample_means ####
    tar_combine(
      seq_accum_sample_means,
      accumulation_plan$accum_map[
        startsWith(names(accumulation_plan$accum_map), "seq_accum_sample_means") &
          !endsWith(names(accumulation_plan$accum_map), "OTU")
      ],
      command = dplyr::bind_rows(!!!.x, .id = "cluster_type") %>%
        dplyr::mutate_at("cluster_type", stringr::str_remove, "seq_accum_sample_means_") %>%
        split_sh_cluster_types()
    ),
    #### sample_accum_plot ####
    tar_target(
      sample_accum_plot,
      # make the plot
      # get the plotting data
      fortify(seq_accum_sample) %>%
        # add forest type annotations
        dplyr::mutate(Type = samples_df[site, "Sites"]) %>%
        # plot it
        ggplot(aes(x = x, y = y, group = cluster_type, color = cluster_type)) +
        # lines for each sample
        geom_line(aes(
          group = paste(site, method, cluster_type),
          # darker when interpolated, lighter when extrapolated
          alpha = dplyr::case_when(method == "interpolated" ~ 0.35, TRUE ~0.10)
        )) +
        # use the alpha values I provided as the literal alpha values
        scale_alpha_identity(name = NULL) +
        scale_color_brewer(type = "qual", palette = 2, name = NULL) +
        # small points to represent the actual observed read depth and richness
        geom_point(data = ~dplyr::filter(., method == "observed"), alpha = 0.4, size = 0.5) +
        # thick dashed line for the mean of all the accumulation curves
        geom_line(data = seq_accum_sample_means, size = 0.8, linetype = "dashed") +
        # seperate facet for each forest type.
        # They are side-by-side so we can compare them visually.
        facet_grid(sh_type~Type, scales = "free_y") +
        xlab("reads") +
        scale_x_continuous() +
        ylab("OTU richness") +
        theme(legend.position = "bottom"),
      packages = c("ggplot2", "iNEXT")
    ),
    #### seq_accum_site ####
    tar_combine(
      seq_accum_site,
      accumulation_plan$accum_map[
        startsWith(names(accumulation_plan$accum_map), "seq_accum_site") &
          !endsWith(names(accumulation_plan$accum_map), "OTU")
      ],
      command = list(!!!.x) %>%
        lapply(ggplot2::fortify) %>%
        dplyr::bind_rows(.id = "cluster_type") %>%
        dplyr::mutate_at("cluster_type", stringr::str_remove, "seq_accum_site_") %>%
        split_sh_cluster_types(),
      packages = "iNEXT"
    ),
    tar_combine(
      samp_accum_site,
      accumulation_plan$accum_map[
        startsWith(names(accumulation_plan$accum_map), "samp_accum_site") &
          !endsWith(names(accumulation_plan$accum_map), "OTU")
      ],
      command = list(!!!.x) %>%
        lapply(ggplot2::fortify) %>%
        dplyr::bind_rows(.id = "cluster_type") %>%
        dplyr::mutate_at("cluster_type", stringr::str_remove, "samp_accum_site_") %>%
        split_sh_cluster_types(),
      packages = "iNEXT"
    ),
    #### asymp_site ####
    tar_combine(
      asymp_site,
      accumulation_plan$accum_map[
        startsWith(names(accumulation_plan$accum_map), "asymp_site") &
          !endsWith(names(accumulation_plan$accum_map), "OTU")
      ],
      command = dplyr::bind_rows(!!!.x, .id = "cluster_type") %>%
        dplyr::mutate_at("cluster_type", stringr::str_remove, "asymp_site_") %>%
        split_sh_cluster_types(),
      packages = "iNEXT"
    ),

    #### *stacked plot* ####
    #### fig_seq_accum_site ####
    tar_target(
      fig_seq_accum_site,
      ggplot(
        seq_accum_site,
        aes(x = x, y = y, ymin = y.lwr, ymax = y.upr,
            # within each facet, each cluster method gets its own line.
            group = cluster_type, color = cluster_type, fill = cluster_type,
            label = cluster_type)
      ) +
        # Dotted horizontal line for the asymptotic estimates
        geom_hline(
          aes(yintercept = Estimator, color = cluster_type),
          linetype = "dashed",
          data = dplyr::filter(asymp_site, xvar == "reads"),
          alpha = 0.5
        ) +
        # ribbon for the confidence band
        geom_ribbon(
          alpha = 0.3, # make it semitransparent
          color = NA # don't draw the edges
        ) +
        # line for the estimated curve
        geom_line(aes(alpha = dplyr::case_when(method == "interpolated" ~ 1, TRUE ~0.5))) +
        scale_alpha_identity(name = NULL) +
        # points for the observed values
        geom_point(data = dplyr::filter(seq_accum_site, method == "observed"),
                   alpha = 1) +
        facet_grid(sh_type ~ site, scales = "free_y") +
        xlab("reads") +
        scale_x_continuous(breaks = 0:3 * 10000,
                           labels = c(0, paste0(1:3 * 10, "k"))) +
        ylab("Species richness") +
        scale_color_brewer(type = "qual", aesthetics = c("color", "fill"),
                           palette = 2, name = NULL) +
        theme(legend.position = "bottom"),
      packages = "ggplot2"
    ),
    #### fig_samp_accum_site ####
    tar_target(
      fig_samp_accum_site,
      ggplot(
        samp_accum_site,
        aes(x = x, y = y, ymin = y.lwr, ymax = y.upr,
            # within each facet, each site gets its own line.
            group = cluster_type, color = cluster_type, fill = cluster_type,
            label = cluster_type)
      ) +
        # Dotted horizontal line for the asymptotic estimates
        geom_hline(
          aes(yintercept = Estimator, color = cluster_type),
          linetype = "dashed",
          data = dplyr::filter(asymp_site, xvar == "samples"),
          alpha = 0.5
        ) +
        # ribbon for the confidence band
        geom_ribbon(
          alpha = 0.3, # make it semitransparent
          color = NA # don't draw the edges
        ) +
        # line for the estimated curve
        geom_line(aes(alpha = dplyr::case_when(method == "interpolated" ~ 1, TRUE ~0.5))) +
        scale_alpha_identity(name = NULL) +
        # points for the observed values
        geom_point(data = dplyr::filter(samp_accum_site, method == "observed"),
                   alpha = 1) +
        facet_grid(sh_type~site, scales = "free_y") +
        xlab("samples") +
        ylab("Species richness") +
        # color brewer has pallettes which are easier to distinguish, although it
        # will always be hard with 8
        scale_color_brewer(type = "qual", aesthetics = c("color", "fill"),
                           palette = 2, name = NULL) +
        theme(legend.position = "bottom"),
      packages = "ggplot2"
    ),
    #### fig_accum_site ####
    tar_target(
      fig_accum_site,
      ggpubr::ggarrange(fig_seq_accum_site, fig_samp_accum_site, ncol = 2,
                        labels = "AUTO", common.legend = TRUE,
                        legend = "bottom")
    ),
    #### *begin map {ext}* ####
    tar_map(
      values = plot_type_meta,
      #### accumplot1_{ext} ####
      tar_file(
        accumplot1,
        write_and_return_file(
          sample_accum_plot,
          file.path(figdir, sprintf("accum1.%s", ext)),
          device = fun, width = 6.25, height = 10, dpi = 150
        )
      ),
      #### accumplot2_{ext} ####
      tar_file(
        accumplot2,
        write_and_return_file(
          fig_accum_site,
          file.path(figdir, sprintf("accum2.%s", ext)),
          device = fun, width = 6.25, height = 5, dpi = 150
        )
      ),
      #### accumplot2_seq_{ext} ####
      tar_file(
        accumplot2_seq,
        write_and_return_file(
          fig_seq_accum_site,
          file.path(figdir, sprintf("accum2_seq.%s", ext)),
          device = fun, width = 4, height = 5, dpi = 150
        )
      ),
      #### accumplot2_samp_{ext} ####
      tar_file(
        accumplot2_samp,
        write_and_return_file(
          fig_samp_accum_site,
          file.path(figdir, sprintf("accum2_samp.%s", ext)),
          device = fun, width = 4, height = 5, dpi = 150
        )
      ),
      names = ext
    ) #### *end map {ext}* ####
  )
)


# # maximum number of reads per sample
# max(colSums(community_matrix))
# min(colSums(community_matrix))
#
#
#
#
#
#
#
#
#
# #What is the greatest number of reads at a condition?
# merge_samples(physeq, "Condition") %>%
#   otu_table() %>%
#   rowSums() %>%
#   max()

