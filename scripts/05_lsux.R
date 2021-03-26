# Run LSUx on the OTU representative sequences

pre_positions_targets <-
  tar_plan(

    # get the path for the CM which is truncated at the LR5 primer site
    # (included in LSUx)
    tar_file(
      cm_32S_trunc,
      system.file(
        file.path("extdata", "fungi_32S_LR5.cm"),
        package = "LSUx"
      )
    ),

    #### Ampliseq clusters ####
    # First load the ampliseq results and use LSUx to cut out subregions
    tar_file(ampliseq_file, "processReads/ampliseq/qiime2_ASV_table.tsv"),

    ampliseq =
      readr::read_tsv(
        ampliseq_file,
        col_types = readr::cols(
          .default = readr::col_double(),
          Feature.ID = readr::col_character(),
          Taxon = readr::col_character(),
          sequence = readr::col_character()
        )
      ) %>%
      dplyr::select("Feature.ID", "sequence") %>%
      tibble::deframe()
  )

positions_targets <-  tar_map(
  values = list(seq = rlang::syms(c("ampliseq", "sl_seqs", "vs_seqs")),
                table = rlang::syms(c("ampliseq_table", "sl_table", "vs_table")),
                id = c("as", "sl", "vs")),
  tar_target(
    positions,
    LSUx::lsux(
      seq = seq[table$OTU],
      cm_32S = cm_32S_trunc,
      ITS1 = TRUE,
      cpu = 8,
      # allow 2 Gb ram (per process)
      mxsize = 2048
    )
  ),

  # Just cut out 5.8S, LSU, ITS2, and ITS
  tar_target(
    regions,
    purrr::map2(
      .x = c("5_8S", "LSU1", "ITS2", "ITS1"),
      .y = c("5_8S", "LSU4", "ITS2", "ITS2"),
      tzara::extract_region,
      seq = seq,
      positions = positions
    ) %>%
      purrr::map2(c("5_8S", "LSU", "ITS2", "ITS"),
                  tibble::enframe, name = "seq_id") %>%
      purrr::reduce(dplyr::full_join, by = "seq_id") %>%
      dplyr::mutate(
        `5_8S_hash` = tzara::seqhash(`5_8S`),
        LSU_hash = tzara::seqhash(LSU),
        ITS2_hash = tzara::seqhash(ITS2),
        ITS_hash = tzara::seqhash(ITS)
      ) %>%
      dplyr::mutate_at(
        dplyr::vars(dplyr::ends_with("_hash")),
        tidyr::replace_na,
        replace = "missing"
      )
  ),
  tar_target(
    hash_key,
    dplyr::left_join(
      dplyr::select(regions, seq_id, dplyr::ends_with("hash")),
      tibble::enframe(tzara::seqhash(seq), name = "seq_id", value = "full_hash"),
      by = "seq_id"
    ) %>%
      dplyr::select(-seq_id) %>%
      unique()
  ),
  names = id
)

lsux_targets <- c(pre_positions_targets, positions_targets)
