# cluster all the ITS2 extracted from the OTU representative sequences

threshold_meta <- tibble::tibble(
  threshold = c(90, 97, 99),
  clustname = c("GH90", "SH97", "SH99")
)

#### *Cluster the ITS2 at different thresholds* ####
its2_cluster_plan <- tar_plan(
  #### its2_file ####
  tar_file(
    its2_file,
    write_and_return_file(
      Biostrings::DNAStringSet(allseqs_ITS2),
      file.path(comparedir, "ITS2.fasta.gz")
    )
  ),
  #### its2_precluster_file ####
  tar_file(
    its2_precluster_file,
    unite_precluster(its2_file, file.path(comparedir, "ITS2_precluster"))
  ),
  #### its2_precluster ####
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
  #### its2_precluster_singletons ####
  its2_precluster_singletons = dplyr::filter(its2_precluster, dplyr::n() == 1),
  #### its2_precluster_clusters ####
  tar_group_by(
    its2_precluster_clusters,
    dplyr::filter(its2_precluster, dplyr::n() > 1),
    cluster
  ),

  tar_map(
    values = threshold_meta,
    names = threshold,
    #### its2_cluster_{threshold} ####
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
      #### cluster_fraction_{type}_{threshold} ####
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
      #### cluster_key_{type}_{threshold} ####
      tar_target(
        cluster_key,
        parse_clusters(its2_cluster, its2_precluster_singletons, "ITS2_hash") %>%
          dplyr::left_join(regions, by = "ITS2_hash") %>%
          dplyr::select(OTU = seq_id, ITS2_hash, !!clustname := cluster),
        tidy_eval = FALSE
      ),
      names = type
    )
  )
)

last_cluster_target <- dplyr::last(its2_cluster_plan)

its2_cluster_plan <- c(
  its2_cluster_plan,
  list(
    #### cluster_key ####
    tar_combine(
      cluster_key,
      last_cluster_target[startsWith(names(last_cluster_target), "cluster_key")],
      command = dplyr::bind_rows(!!!.x) %>%
        tidyr::pivot_longer(cols = c("GH90", "SH97", "SH99"), names_to = "SH", values_to = "clust") %>%
        dplyr::filter(complete.cases(.)) %>%
        tidyr::pivot_wider(names_from = "SH", values_from = "clust")
    ),
    #### cluster_detect_table ####
    tar_fst_tbl(
      cluster_detect_table,
      tidyr::pivot_longer(cluster_key, cols = GH90:SH99, names_to = "SH",
                          values_to = "clust") %>%
        dplyr::filter(complete.cases(.)) %>%
        dplyr::mutate(
          SH = as.numeric(substr(SH, 3, 4)),
          OTU = substr(OTU, 1, 3)
        ) %>%
        dplyr::group_by(SH, clust) %>%
        dplyr::summarize(n = dplyr::n_distinct(OTU), .groups = "drop_last") %>%
        dplyr::filter(n == 1) %>%
        dplyr::summarize("Only one method" = dplyr::n(), .groups = "drop")
    ),
    #### cluster_fraction_table ####
    tar_combine(
      cluster_fraction_table,
      last_cluster_target[startsWith(names(last_cluster_target), "cluster_fraction")],
      command = dplyr::bind_rows(!!!.x) %>%
        recode_cluster_types() %>%
        tidyr::pivot_wider(names_from = cluster_type, values_from = fraction) %>%
        dplyr::left_join(cluster_detect_table, by = c("threshold" = "SH"))
    ),
    #### write_cluster_fraction_table ####
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
