# cluster all the ITS2 extracted from the OTU representative sequences

#### Cluster the ITS2 at different thresholds ####
its2_cluster_targets <- tar_plan(
  tar_file(
    its2_file,
    file.path(comparedir, "ITS2.fasta.gz") %T>%
      Biostrings::writeXStringSet(Biostrings::DNAStringSet(allseqs_ITS2), .)
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
    )
  )
)
