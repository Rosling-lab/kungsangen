deseq_plan <-
  tar_plan(
    ########################################################################
    ########################################################################
    #otu_a
    phyloseq_otua =
      phyloseq::merge_phyloseq(
        samples_physeq_fungi_otuA,
        samples_physeq_protists_otuA
      ) %>%
      phyloseq::subset_samples(Sites != "None"),
    cds_otua = phyloseq::phyloseq_to_deseq2(phyloseq_otua, ~Sites) %>%
      DESeq2::DESeq(test = "Wald", fitType = "parametric"),

    res_otua = DESeq2::results(cds_otua, cooksCutoff = TRUE),
    fullres_otua = DESeq2::results(cds_otua, cooksCutoff = FALSE),

    alpha = 0.05,
    sigtab_otua = as(res_otua, "data.frame") %>% {
      .[which(.$padj < alpha),]
    } %>% {
      cbind(
        as(., "data.frame"),
        as(phyloseq::tax_table(phyloseq_otua)[rownames(.),], "matrix")
      )
    },
    tar_file(
      name = sigtab_otua_file,
      file.path(datadir, "sigtab_otua.txt") %T>%
        write.table(sigtab_otua, ., row.names =T , col.names = T, sep = "\t")
    ),
    # function defined in 01_functions.R
    deseq_plot_otua = deseq_plot(sigtab_otua, "OTU_A"),
    # save figure
    tar_map(
      values = plot_type_meta,
      names = ext,
      tar_file(
        name = deseq_plot_otua_file,
        write_and_return_file(
          deseq_plot_otua,
          file.path(figdir, paste0("deseq_plot_otua.", ext)),
          device = fun, width = 4.5, height = 4, dpi = 150
        )
      )
    ),

    ########################################################################
    ########################################################################
    #otu_s
    phyloseq_otus =
      phyloseq::merge_phyloseq(
        samples_physeq_fungi_otuS,
        samples_physeq_protists_otuS
      ) %>%
      phyloseq::subset_samples(Sites != "None"),
    cds_otus = phyloseq::phyloseq_to_deseq2(phyloseq_otus,  ~Sites) %>%
      DESeq2::DESeq(test = "Wald", fitType = "parametric"),

    res_otus= DESeq2::results(cds_otus, cooksCutoff = TRUE),

    fullres_otus = DESeq2::results(cds_otus, cooksCutoff = FALSE),
    sigtab_otus = as(res_otus, "data.frame") %>% {
      .[which(.$padj < alpha),]
    } %>% {
      cbind(
        as(., "data.frame"),
        as(phyloseq::tax_table(phyloseq_otus)[rownames(.), ], "matrix")
      )
    },
    tar_file(
      sigtab_otus_file,
      file.path(datadir, "sigtab_otus.txt") %T>%
        write.table(sigtab_otus, ., row.names =T , col.names = T, sep = "\t")
    ),
    deseq_plot_otus = deseq_plot(sigtab_otus, "OTU_S"),
    # save figure
    tar_map(
      values = plot_type_meta,
      names = ext,
      tar_file(
        deseq_plot_otus_file,
        write_and_return_file(
          deseq_plot_otus,
          file.path(figdir, paste0("deseq_plot_otus.", ext)),
          device = fun, width = 4.5, height = 4, dpi = 150
        )
      )
    ),

    #otu_c
    phyloseq_otuc =
      phyloseq::merge_phyloseq(
        samples_physeq_fungi_otuC,
        samples_physeq_protists_otuC
      ) %>%
      phyloseq::subset_samples(Sites != "None"),
    cds_otuc = phyloseq::phyloseq_to_deseq2(phyloseq_otuc,  ~Sites) %>%
      DESeq2::DESeq(test = "Wald", fitType = "parametric"),

    res_otuc = DESeq2::results(cds_otuc, cooksCutoff = TRUE),

    fullres_otuc = DESeq2::results(cds_otuc, cooksCutoff = FALSE),
    sigtab_otuc= as(res_otuc, "data.frame") %>% {
      .[which(.$padj < alpha),]
    } %>% {
      cbind(
        as(., "data.frame"),
        as(phyloseq::tax_table(phyloseq_otuc)[rownames(.), ], "matrix")
      )
    },
    tar_file(
      sigtab_otuc_file,
      file.path(datadir, "sigtab_otuc.txt") %T>%
        write.table(sigtab_otuc, ., row.names =T , col.names = T, sep = "\t")
    ),
    # function defined in 01_functions.R
    deseq_plot_otuc = deseq_plot(sigtab_otuc, "OTU_C"),
    # save figure
    tar_map(
      values = plot_type_meta,
      names = ext,
      tar_file(
        deseq_plot_otuc_file,
        write_and_return_file(
          deseq_plot_otuc,
          file.path(figdir, paste0("deseq_plot_otuc.", ext)),
          device = fun, width = 4.5, height = 4, dpi = 150
        )
      )
    )
  )
