# save OTU tables for just fungi and just "protists"

# map across taxonomic group, cluster type, and file format
tar_map(
  values = list(
    group = c("protists", "fungi"),
    tree = rlang::syms(c("tree_fungi", "tree_protists"))
  ),
  tar_map(
    values = list(
      clust_type = c("ampliseq", "vsearch", "swarm"),
      otutab = rlang::syms(c("ampliseq_table", "vs_table", "sl_table")),
      regions = rlang::syms(c("regions_as", "regions_vs", "regions_sl"))
    ),
    tar_fst_tbl(
      otu_table,
      dplyr::select(regions, OTU = seq_id, `5_8S_hash`, LSU_hash) %>%
        tidyr::unite("label", `5_8S_hash`, LSU_hash) %>%
        dplyr::filter(label %in% tree$tip.label) %>%
        dplyr::left_join(otutab, by = "OTU") %>%
        dplyr::select(-label)
    ),
    tar_map(
      values = list(
        ext = c("rds", "xlsx")
      ),
      tar_file(
        write_otu_table,
        file.path(datadir, sprintf("%s_table_%s.%s", clust_type, group, ext)) %>%
          write_table(otu_table, file = ., format = ext)
      )
    ),
    names = clust_type
  ),
  names = group
) -> kingdom_table_targets
