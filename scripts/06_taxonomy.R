# Assign taxonomy

cumcollapse <- function(x) {
  purrr::map_chr(
    seq_along(x),
    ~paste(x[1:.], collapse = ";")
  )
}

## identifies a clade which is the outgroup based on the confident
# kingdom assignments, and returns all the tips in it
root_with_kingdoms <- function(tree, kingdoms, outgroup_list) {
  outgroup <- dplyr::filter(kingdoms, taxon %in% outgroup_list)
  ingroup <- dplyr::filter(kingdoms, !taxon %in% outgroup_list)
  outgroup_mrca <- ape::getMRCA(tree, intersect(outgroup$label, tree$tip.label))
  # if the root is currently inside outgroup, root it with fungi first instead
  if (outgroup_mrca == phangorn::getRoot(tree)) {
    fungi <- dplyr::filter(kingdoms, taxon %in% "Fungi")
    nonfungi <- dplyr::filter(kingdoms, !taxon %in% "Fungi")
    fungi_mrca <- ape::getMRCA(tree, intersect(fungi$label, tree$tip.label))
    if (fungi_mrca == phangorn::getRoot(tree)) {
      stop("Fungi and outgroup (", outgroup_list,
           ") are paraphyletic with respect to each other, ",
           "cannot root tree.")
    }
    fungi_tips <- unlist(phangorn::Descendants(tree, fungi_mrca))
    tree <- ape::root.phylo(tree, fungi_tips, edgelabel = TRUE)
    outgroup_mrca <- ape::getMRCA(tree, outgroup$label)
  }
  outgroup <- tree$tip.label[unlist(phangorn::Descendants(tree, outgroup_mrca))]
  # error if outgroup is not monophyletic with respect to the other
  # confident assignments
  stopifnot(!any(ingroup$label %in% outgroup))
  ape::root.phylo(tree, outgroup, resolve.root = TRUE, edgelabel = TRUE)
}

extract_fungi <- function(kingdoms, tree) {
  fungi <- dplyr::filter(kingdoms, taxon == "Fungi")
  nonfungi <- dplyr::filter(kingdoms, taxon != "Fungi")
  fungi_mrca <- ape::getMRCA(tree, fungi$label)
  nonfungi_mrca <- ape::getMRCA(tree, nonfungi$label)
  if (fungi_mrca == phangorn::getRoot(tree)) {
    # the tree is currently rooted inside fungi, so re-root it using all non-fungi
    ape::root.phylo(
      tree,
      tree$tip.label[]
    )
  }
}

check_monophyly <- function(taxtable, tree, taxa) {
  tips <- dplyr::filter(taxtable, taxon %in% taxa)$label
  mrca <- ape::getMRCA(tree, tips)
  clade_tips <- tree$tip.label[unlist(phangorn::Descendants(tree, mrca))]
  bad_taxa <- dplyr::filter(taxtable, label %in% clade_tips, !taxon %in% taxa)$taxon
  if (length(bad_taxa) > 0) {
    stop("Clade (", paste(taxa, collapse = ", "), ") is polyphyletic with respect to:\n",
         glue::glue("-{names(tab)}: {tab} sequences\n", tab = table(bad_taxa)))
  }
}

check_tree <- function(kingdoms, tree) {
  lapply(
    c(
      as.list(unique(kingdoms$taxon)),
      list(
        c("Stramenopila", "Alveolata", "Rhizaria"),
        c("Stramenopila", "Alveolata", "Rhizaria", "Viridiplantae"),
        c("Fungi", "Metazoa"),
        c("Fungi", "Metazoa", "Amoebozoa")
      )
    ),
    check_monophyly,
    taxtable = kingdoms,
    tree = tree
  )
}

# take a subtree defined by a clade name and a taxonomy table
extract_clade <- function(tree, clade, taxtable) {
  tips <- dplyr::filter(taxtable, taxon %in% clade)
  mrca <- ape::getMRCA(tree, tips)
  clade_tips <- unlist(phangorn::Descendants(tree, mrca))
  ape::keep.tip(tree, clade_tips)
}

taxon_abbrevs <- tibble::tribble(
  ~pattern,      ~replacement,
  "(mycota|mycetes|ales|aceae|plantae|zoa|phyceae|phyta|phora|monada)(\\b|_)", "\\2",
  "Fungi(\\b|_)",    "Fun\\1",
  "Basidio(\\b|_)",  "Bas\\1",
  "Asco(\\b|_)",     "Asc\\1",
  "Chytridio(\\b|_)","Chy\\1",
  "Mucoro(\\b|_)",   "Muc\\1",
  "Mortierello(\\b|_)","Mor\\1",
  "Zoopago(\\b|_)",  "Zpg\\1",
  "Kickxello(\\b|_)","Kix\\1",
  "Monoblepharo(\\b|_)","Mnb\\1",
  "Entorrhizo(\\b|_)","Etr\\1",
  "Entomophthoro(\\b|_)","Etm\\1",
  "Calcarisporiello(\\b|_)","Cal\\1",
  "Blastocladio(\\b|_)","Bla\\1",
  "Neocallimastigo(\\b|_)","Neo\\1",
  "Glomero(\\b|_)",  "Glo\\1",
  "Olpidio(\\b|_)",  "Olp\\1",
  "Rhizaria(\\b|_)", "Rhi\\1",
  "Stramenopila(\\b|_)","Stram\\1",
  "Alveolata(\\b|_)", "Alv\\1",
  "incertae_sedis", "?",
  "unidentified", "??"
)

taxonomy_targets <- tar_plan(
  tar_map(
    values = list(
      allseqs = rlang::syms(c("allseqs_ITS", "allseqs_LSU", "allseqs_LSU")),
      reference = c("unite.ITS", "rdp_train.LSU", "silva.LSU"),
      id = c("unite", "rdp", "silva")
    ),
    tar_file(
      reference_file,
      sprintf("reference/%s.sintax2.fasta.gz", reference)
    ),
    # unique kingdoms in each reference dataset
    tar_qs(
      ref_kingdoms,
      system(
        paste("zcat", reference_file, "| sed -nr 's/.*k:([^,]*).*/\\1/p' | sort -u"),
        intern = TRUE
      )
    ),
    tar_fst_tbl(
      tax,
      phylotax::taxonomy_sintax(
        seq = allseqs,
        reference = reference_file,
        min_confidence = 80,
        multithread = 8
      ) %>%
        phylotax::taxtable_sintax()
    ),
    names = id
  ),
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
  tar_file(
    constraint_file,
    file.path("reference", "constraints.fasta")
  ),
  tax_constraints = Biostrings::readBStringSet(constraint_file) %>%
    as.character() %>%
    tibble::enframe(name = "taxon", value = "constraint") %>%
    dplyr::left_join(dplyr::select(kingdoms, taxon, label), by = "taxon") %>%
    dplyr::select(label, constraint) %>%
    dplyr::filter(complete.cases(.)) %>%
    tibble::deframe() %>%
    Biostrings::BStringSet()
)

taxonomy_targets <- c(
  taxonomy_targets,
  tar_combine(
    kingdom_count,
    taxonomy_targets[[1]]$ref_kingdoms,
    command = table(c(!!!.x))
  )
)
