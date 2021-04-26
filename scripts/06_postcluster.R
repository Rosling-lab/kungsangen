# cluster all the ITS2 extracted from the OTU representative sequences

#### Cluster the ITS2 at different thresholds ####
its2_cluster_targets <- tar_plan(
  tar_file(
    its2_file,
    write_and_return_file(
      Biostrings::DNAStringSet(allseqs_ITS2),
      file.path(comparedir, "ITS2.fasta.gz")
    )
  ),

  tar_file(
    its2_precluster_file,
    unite_precluster(its2_file, file.path(comparedir, "ITS2_precluster"))
  ),
  its2_precluster = readr::read_tsv(
    its2_precluster_file,
    col_names = paste0("V", 1:10),
    col_types = "ciidfccccc",
    na = c("", "NA", "*")
  ) %>%
    dplyr::select(seq_id = V9, cluster = V10) %>%
    dplyr::mutate_all(stringr::str_replace, ";.*", "") %>%
    dplyr::mutate(cluster = dplyr::coalesce(cluster, seq_id)) %>%
    unique() %>%
    dplyr::group_by(cluster),
  its2_precluster_singletons = dplyr::filter(its2_precluster, dplyr::n() == 1),
  tar_group_by(
    its2_precluster_clusters,
    dplyr::filter(its2_precluster, dplyr::n() > 1),
    cluster
  ),

  tar_map(
    values = list(threshold = c(90, 97, 99)),
    tar_target(
      its2_cluster,
      unite_cluster(its2_precluster_clusters, allseqs_ITS2, threshold),
      pattern = map(its2_precluster_clusters)
    ),
    tar_map(
      values = list(
        regions = rlang::syms(c("regions_as", "regions_vs", "regions_sl")),
        type = c("as", "vs", "sl")
      ),
      tar_target(
        cluster_fraction,
        parse_clusters(its2_cluster, its2_precluster_singletons, "ITS2_hash") %>%
          dplyr::left_join(regions, by = "ITS2_hash") %>%
          dplyr::group_by(cluster) %>%
          dplyr::summarize(present = !all(is.na(seq_id))) %>%
          dplyr::summarize(
            cluster_type = type,
            threshold = threshold,
            total = dplyr::n(),
            fraction = sum(present)/dplyr::n()
          )
      ),
      names = type
    )
  )
)

last_cluster_target <- dplyr::last(its2_cluster_targets)

its2_cluster_targets <- c(
  its2_cluster_targets,
  list(
    tar_combine(
      cluster_fraction_table,
      last_cluster_target[startsWith(names(last_cluster_target), "cluster_fraction")],
      command = dplyr::bind_rows(!!!.x) %>%
        recode_cluster_types() %>%
        tidyr::pivot_wider(names_from = cluster_type, values_from = fraction)
    ),
    tar_file(
      write_cluster_fraction_table,
      write_and_return_file(
        cluster_fraction_table,
        file.path(datadir, "cluster_fractions.xlsx"),
        "xlsx"
      )
    )
  )
)
