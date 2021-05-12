# Run LSUx on the OTU representative sequences

# get a linked copy of the truncated 32S CM from LSUx
lsux_cm_32S <- system.file(
  file.path("extdata", "fungi_32S_LR5.cm"),
  package = "LSUx"
)
cm_32S_trunc_file <- file.path("reference", "fungi_32S_LR5.cm")
if (file.exists(cm_32S_trunc_file) &&
    tools::md5sum(cm_32S_trunc_file) != tools::md5sum(lsux_cm_32S)) {
  unlink(cm_32S_trunc_file)
}
if (!file.exists(cm_32S_trunc_file)) {
  file.link(
    lsux_cm_32S,
    cm_32S_trunc_file
  )
}

pre_positions_plan <-
  tar_plan(

    # get the path for the CM which is truncated at the LR5 primer site
    # (included in LSUx)
    tar_file(
      cm_32S_trunc,
      cm_32S_trunc_file
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

positions_meta <- tibble::tibble(
  seq = rlang::syms(c("ampliseq", "seqs_sl", "seqs_vs")),
  table = rlang::syms(c("ampliseq_table", "table_sl", "table_vs")),
  id = c("as", "sl", "vs")
)

positions_plan <-  tar_map(
  values = positions_meta,
  tar_target(
    positions,
    LSUx::lsux(
      seq = seq[table$OTU],
      cm_32S = cm_32S_trunc,
      ITS1 = TRUE,
      cpu = local_cpus(),
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
      dplyr::full_join(
        tibble::enframe(seq, name = "seq_id", value = "full"),
        by = "seq_id"
      ) %>%
      dplyr::mutate(
        `5_8S_hash` = tzara::seqhash(`5_8S`),
        LSU_hash = tzara::seqhash(LSU),
        ITS2_hash = tzara::seqhash(ITS2),
        ITS_hash = tzara::seqhash(ITS),
        full_hash = tzara::seqhash(full)
      ) %>%
      dplyr::mutate_at(
        dplyr::vars(dplyr::ends_with("_hash")),
        tidyr::replace_na,
        replace = "missing"
      ) %>%
      tidyr::unite("label", "5_8S_hash", "LSU_hash", remove = FALSE)
  ),
  tar_target(
    hash_key,
    dplyr::select(regions, seq_id, label, dplyr::ends_with("hash")) %>%
      dplyr::rename(!!id := seq_id) %>%
      unique(),
    tidy_eval = FALSE
  ),
  names = id
)

lsux_plan <- c(pre_positions_plan, positions_plan)
