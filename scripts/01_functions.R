# Functions to use in targets pipeline

#### Functions for loading data ####

load_cons_seqs <- function(seq_file, prefix) {
  seqs <- Biostrings::readDNAStringSet(seq_file)
  names(seqs) <- gsub("consensus", prefix, names(seqs))
  as.character(seqs)
}

load_sl_rawtable <- function(sl_table_file, sl_seqs) {
  readr::read_delim(
    sl_table_file,
    delim = " ",
    col_names = c("Sample", "OTU", "reads"),
    col_types = "cci"
  ) %>%
    dplyr::arrange(OTU) %>%
    tidyr::pivot_wider(
      names_from = "Sample",
      values_from = "reads",
      values_fill = 0L
    ) %>%
    dplyr::left_join(
      tibble::enframe(sl_seqs, name = "OTU", value = "seq"),
      by = "OTU"
    )
}

#### Functions for Taxonomy assignment ####
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

#### Tree building ####

# make an (optionally constrainted) ML tree with fasttree
fasttree <- function(aln_file, out_file, constraints = NULL) {
  args <- character(0)
  if (!is.null(constraints)) {
    constraint_file <- tempfile(fileext = ".fasta")
    Biostrings::writeXStringSet(constraints, constraint_file)
    on.exit(unlink(constraint_file))
    args <- c("-constraints", constraint_file)
  }
  stopifnot(
    system2(
      command = "fasttree",
      args = c(args, "-nt", "-gtr", aln_file),
      stdout = out_file
    ) == 0
  )
  out_file
}

#### Unite SH clustering ####

unite_precluster <- function(seqs, outroot) {
  assertthat::assert_that(
    assertthat::is.string(seqs),
    assertthat::is.readable(seqs)
  )
  out_uc = glue::glue("{outroot}.uc")
  stopifnot(
    system2(
      command = "vsearch",
      args = c("--cluster_fast", seqs, "--id", "0.80", "--uc", out_uc)
    ) == 0
  )
  out_uc
}

unite_cluster <- function(preclust, seqs, threshold, outfile) {
  seqs <- Biostrings::DNAStringSet(seqs[preclust$seq_id])
  tf <- tempfile(pattern = "clust", fileext = ".fasta")
  outfile <- tempfile(pattern = "out")
  Biostrings::writeXStringSet(seqs, tf)
  on.exit(unlink(tf))
  on.exit(unlink(outfile), TRUE)
  stopifnot(
    system2(
      "blastclust",
      c("-i", tf, "-S", threshold, "-L", "0.85", "-a", "8", "-e", "F", "-o",
        outfile, "-p", "F")
    ) == 0
  )
  readLines(outfile)
}

#### Displaying and annotating trees ####


draw_cluster <- function(node, offset, extend = 0.35, barsize = 2, ...) {
  geom_cladelab(node, "", align = TRUE, offset = offset, extend = extend, barsize = barsize, ...)
}

draw_clusters <- function(clusters, singletons, hash_key, physeq, offset, name = "") {
  tree <- phyloseq::phy_tree(physeq)
  tips <- phyloseq::taxa_names(physeq)
  clusters <- dplyr::bind_rows(
    tibble::tibble(
      seq_id = trimws(clusters),
      cluster = formatC(seq_along(clusters), width = 4, flag = "0")
    ),
    singletons
  ) %>%
    dplyr::rename(ITS2_hash = seq_id) %>%
    tidyr::separate_rows(ITS2_hash) %>%
    dplyr::left_join(hash_key, by = "ITS2_hash") %>%
    dplyr::transmute(cluster = cluster, label = paste(`5_8S_hash`, LSU_hash, sep = "_")) %>%
    dplyr::filter(label %in% tips) %>%
    unique() %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarize(
      n = dplyr::n(),
      mrca = if (n > 1) ape::getMRCA(tree, label) else  match(label, tips),
      n_desc = length(unlist(phangorn::Descendants(tree, mrca[1])))
    )
  mono_clusters <- dplyr::filter(clusters, n == n_desc)
  poly_clusters <- dplyr::filter(clusters, n < n_desc)
  list(
    draw_cluster(poly_clusters$mrca, offset = offset, barcolor = scales::alpha("red", 0.5)),
    draw_cluster(mono_clusters$mrca, offset = offset, barsize = 1.5)
  )
}
cluster_annotation <- function(label, x, y = -0.20) {
  ggplot2::annotation_custom(
    grid::textGrob(label = label, rot = 90, hjust = 1, gp = grid::gpar(fontsize = 6)),
    xmin = x, xmax = x, ymin = y, ymax = y)
}

#### Functions for species accumulation curves ####

# function to print numbers greater than 1000 with "k"
formatk <- function(x, ...) {
  ifelse(
    x < 1000,
    format(x, ...),
    paste0(format(x / 1000, ...), "k")
  )
}

# function to make mean of species accumulation curves
# accum is the result of iNEXT::iNEXT
# new_x is the points to estimate the mean of the accumulation curves
# samples_data is a data.frame or matrix with the same rownames used in the
#   otu table passed to iNEXT
# groupname is the name of a column in samples_data which defines groups of sites
#   to average over
accumulation_means <- function(accum, newx, sampledata, groupname) {
# get the plotting data
  ggplot2::fortify(accum) %>%
    # iNEXT doesn't give the same set of x values for each site, so it's hard to
    # take a reasonable average.  So for each site, make a spline function to
    # interpolate the iNEXT results,and get the values at a particular set of x
    # values.
    dplyr::group_by(site) %>%
    dplyr::summarize(
      interp = list(splinefun(x = x, y = y)),
      x = list(newx),
      y = list(interp[[1]](x[[1]]))
    ) %>%
    dplyr::select(site, x, y) %>%
    tidyr::unchop(c(x, y)) %>%
    # add wetness annotations
    dplyr::mutate(Type = sampledata[as.character(site), groupname]) %>%
    # get the median at each wetness.
    dplyr::group_by(x, Type) %>%
  dplyr::summarize_at("y", mean, .groups = "none")
}

recode_cluster_types <- function(data) {
  dplyr::mutate_at(
    data,
    "cluster_type",
    dplyr::recode,
    as = "ASV (Ampliseq/DADA2)",
    sl = "Single-linkage (GeFaST)",
    vs = "Centroid-based (VSEARCH)"
  )
}


# Save a data frame as XLSX or RDS
write_table <- function(table, file, format = c("rds", "xlsx")) {
  format <- match.arg(format)
  switch(
    format,
    rds = saveRDS(object = table, file = file, version = 2),
    xlsx = openxlsx::write.xlsx(x = table, file = file)
  )
  file
}
