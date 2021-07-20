# save OTU tables for just fungi and just "protists"

group_meta <- tibble::tibble(
  group = c("protists", "fungi"),
  tree = rlang::syms(c("tree_protists", "tree_fungi_new"))
)

table_meta <- tibble::tibble(
  clust_type = c("ampliseq", "vsearch", "swarm"),
  otutab = rlang::syms(c("ampliseq_table", "table_vs", "table_sl")),
  regions = rlang::syms(c("regions_as", "regions_vs", "regions_sl"))
)

kingdom_table_meta <-
  dplyr::full_join(group_meta, table_meta, by = character()) %>%
  dplyr::mutate(
    kingtab = rlang::syms(sprintf("otu_table_%s_%s", clust_type, group))
  )

# map across taxonomic group, cluster type, and file format
list(
  tar_map(
    values = group_meta,
    names = group,
    tar_map(
      values = table_meta,
      names = clust_type,
      #### otu_table_{clust_type}_{group} ####
      tar_fst_tbl(
        otu_table,
        dplyr::select(regions, OTU = seq_id, label) %>%
          dplyr::filter(label %in% tree$tip.label) %>%
          dplyr::left_join(otutab, by = "OTU") %>%
          dplyr::select(-label)
      ),
      tar_map(
        values = list(
          ext = c("rds", "xlsx")
        ),
        #### write_otu_table_{ext}_{clust_type}_{group} ####
        tar_file(
          write_otu_table,
          file.path(datadir, sprintf("%s_table_%s.%s", clust_type, group, ext)) %>%
            write_table(otu_table, file = ., format = ext)
        )
      )
    )
  ),
  tar_target(
    kingdom_table,
    tibble::tibble(
      clust_type = !!kingdom_table_meta$clust_type,
      group = !!kingdom_table_meta$group,
      otutab = list(!!!kingdom_table_meta$otutab),
      kingtab = list(!!!kingdom_table_meta$kingtab)
    ) %>%
      dplyr::mutate_at(
        c("otutab", "kingtab"),
        purrr::map,
        tibble::column_to_rownames,
        "OTU"
      ) %>%
      dplyr::group_by(clust_type) %>%
      dplyr::transmute(
        clust_type = clust_type,
        group = group,
        tot_reads = purrr::map_int(otutab, sum),
        reads = purrr::map_int(kingtab, sum),
        micro_reads = sum(reads),
        tot_otus = purrr::map_int(otutab, nrow),
        otus = purrr::map_int(kingtab, nrow),
        micro_otus = sum(otus),
        frac_reads = reads/tot_reads,
        frac_otus = otus/tot_otus,
        frac_micro_reads = reads/micro_reads,
        frac_micro_otus = otus/micro_otus
      ),
    tidy_eval = TRUE
  )
) -> kingdom_table_plan
