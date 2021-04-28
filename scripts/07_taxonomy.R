# Assign taxonomy

taxonomy_targets <- tar_plan(
  tar_map(
    values = list(
      # allseqs is actually defined in 08_tree.R, but that's ok, targets will
      # execute them in the correct order.
      allseqs = rlang::syms(c("allseqs_ITS", "allseqs_LSU", "allseqs_LSU")),
      reference = c("unite_single.ITS", "rdp_train.LSU", "silva_nr99.LSU"),
      id = c("unite", "rdp", "silva")
    ),
    # files for the reference datasets
    tar_file(
      reference_file,
      sprintf("reference/%s.sintax.fasta.gz", reference)
    ),
    # unique kingdoms in each reference dataset
    tar_qs(
      ref_kingdoms,
      system(
        paste("zcat", reference_file, "| sed -nr 's/.*k:([^,]*).*/\\1/p' | sort -u"),
        intern = TRUE
      )
    ),
    # assign taxonomy using sintax
    tar_fst_tbl(
      tax,
      phylotax::taxonomy_sintax(
        seq = allseqs,
        reference = reference_file,
        min_confidence = 80,
        multithread = local_cpus()
      ) %>%
        phylotax::taxtable_sintax()
    ),
    names = id
  ),
  # combine all the taxonomy results into one table
  tax_all = dplyr::bind_rows(
    tibble::add_column(tax_unite, method = "unite") %>%
      dplyr::rename(ITS_hash = label) %>%
      dplyr::left_join(dplyr::select(hash_key, ITS_hash, full_hash),
                       by = "ITS_hash"),
    tibble::add_column(tax_rdp, method = "rdp") %>%
      dplyr::rename(LSU_hash = label) %>%
      dplyr::left_join(dplyr::select(hash_key, LSU_hash, full_hash),
                       by = "LSU_hash"),
    tibble::add_column(tax_silva, method = "silva") %>%
      dplyr::rename(LSU_hash = label) %>%
      dplyr::left_join(dplyr::select(hash_key, LSU_hash, full_hash),
                       by = "LSU_hash")
  ) %>%
    dplyr::filter(confidence >= 0.8),
  # for each taxonomic level, when there are multiple assignments, do they match?
  tax_both =
    dplyr::filter(
        tax_all,
        !startsWith(taxon, "unidentified"),
        !endsWith(taxon, "ncertae_sedis")
    ) %>%
    dplyr::group_by(full_hash, rank) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::summarize(match = dplyr::n_distinct(taxon) == 1) %>%
    dplyr::ungroup() %>%
    dplyr::add_count(rank) %>%
    dplyr::filter(match) %>%
    dplyr::group_by(rank, n) %>%
    dplyr::summarize(frac = dplyr::n() / unique(n)),
  # what are the disagreements between the different methods?
  tax_compare =
    tax_all %>%
    dplyr::group_by(full_hash, method) %>%
    dplyr::arrange(rank) %>%
    dplyr::mutate(c12n = cumcollapse(taxon)) %>%
    dplyr::group_by(full_hash, rank) %>%
    dplyr::filter(
      !startsWith(taxon, "unidentified"),
      !endsWith(taxon, "ncertae_sedis")
    ) %>%
    dplyr::filter(dplyr::n_distinct(taxon) > 1) %>%
    dplyr::select(full_hash, taxon, rank, method, c12n) %>%
    tidyr::pivot_wider(names_from = "method", values_from = c("taxon", "c12n")) %>%
    dplyr::ungroup() %>%
    dplyr::count(rank, taxon_unite, taxon_rdp, taxon_silva, c12n_unite,
                  c12n_rdp, c12n_silva),

  # get the kingdom assignments where there are no conflicts
  # these will be used to constrain the tree.
  kingdoms = tax_all %>%
    # take kingdom assignments
    dplyr::filter(rank == "kingdom") %>%
    dplyr::group_by(full_hash) %>%
    # keep only assignments where there is no disagreement
    dplyr::filter(dplyr::n_distinct(taxon) == 1) %>%
    dplyr::select(taxon, full_hash) %>%
    unique() %>%
    # join to the has key to get hashes of 5.8S and LSU
    # since these are the tip names on the tree
    dplyr::left_join(hash_key, by = "full_hash") %>%
    dplyr::mutate(label = paste(`5_8S_hash`, LSU_hash, sep = "_")) %>%
    dplyr::group_by(label) %>%
    # require unanimity among identifications
    dplyr::filter(dplyr::n_distinct(taxon) == 1),
  # constraints for the kingdoms
  # enforces monophyly of each kingdom, as well as supergroups
  tar_file(
    constraint_file,
    file.path("reference", "constraints.fasta")
  ),
  # copy the constraints for each kingdom onto the sequences which
  # belong to that kingdom
  tax_constraints = Biostrings::readBStringSet(constraint_file) %>%
    as.character() %>%
    tibble::enframe(name = "taxon", value = "constraint") %>%
    dplyr::left_join(dplyr::select(kingdoms, taxon, label), by = "taxon") %>%
    dplyr::select(label, constraint) %>%
    dplyr::filter(complete.cases(.)) %>%
    tibble::deframe() %>%
    Biostrings::BStringSet(),
  # get the plylum assignments where there are no conflicts
  # these are only used to pick out the vascular plants
  phyla = tax_all %>%
    # take kingdom and phylum assignments
    dplyr::filter(rank %in% c("kingdom", "phylum")) %>%
    # deal separately with kindoms and phyla
    dplyr::group_split(rank) %>%
    purrr::map(dplyr::group_by, full_hash) %>%
    # keep only assignments where there is no disagreement
    purrr::map(dplyr::filter, dplyr::n_distinct(taxon) == 1) %>%
    {dplyr::semi_join(.[[2]], .[[1]], by = "full_hash")} %>%
    dplyr::select(taxon, full_hash) %>%
    unique() %>%
    # join to the has key to get hashes of 5.8S and LSU
    # since these are the tip names on the tree
    dplyr::left_join(hash_key, by = "full_hash") %>%
    dplyr::mutate(label = paste(`5_8S_hash`, LSU_hash, sep = "_")) %>%
    dplyr::group_by(label) %>%
    # require unanimity among identifications
    dplyr::filter(dplyr::n_distinct(taxon) == 1)
)

taxonomy_targets <- c(
  taxonomy_targets,
  tar_combine(
    kingdom_count,
    taxonomy_targets[[1]]$ref_kingdoms,
    command = table(c(!!!.x))
  )
)
