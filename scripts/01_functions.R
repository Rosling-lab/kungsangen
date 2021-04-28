# Functions to use in targets pipeline

#### Utility functions ####

# are we running slurm?
is_slurm <- function() nchar(Sys.which("sbatch")) > 0
is_local <- function() !is_slurm()

# are we running snakemake?
is_snakemake <- function() !interactive() && exists("snakemake")

# how many cpus do we have on the local machine?
# if we're not running on the cluster, leave one cpu free.
local_cpus <- function() {
  if (is_snakemake()) {
    snakemake@threads
  } else if (is_slurm()) {
    out <- as.integer(Sys.getenv("SLURM_JOB_CPUS_PER_NODE"))
    assertthat::assert_that(assertthat::is.count(out))
    out
  } else {
    max(parallel::detectCores() - 1, 1)
  }
}

#### Functions for loading data ####

load_cons_seqs <- function(seq_file, prefix) {
  seqs <- Biostrings::readDNAStringSet(seq_file)
  names(seqs) <- gsub("consensus", prefix, names(seqs))
  as.character(seqs)
}

load_rawtable <- function(table_file, seqs) {
  readr::read_delim(
    table_file,
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
      tibble::enframe(seqs, name = "OTU", value = "seq"),
      by = "OTU"
    )
}

#### convenience function for writing a file and returning its name ####

ensure_directory <- function(file) {
  d <- dirname(file)
  if (!dir.exists(d)) dir.create(d)
}

write_and_return_file <- function(x, file, ...) {
  UseMethod("write_and_return_file")
}

write_and_return_file.XStringSet <- function(x, file, ...) {
  ensure_directory(file)
  Biostrings::writeXStringSet(x, file, ...)
  file
}

write_and_return_file.data.frame <- function(x, file, type = c("rds", "xlsx"), ...) {
  ensure_directory(file)
  type = match.arg(type)
  switch(
    type,
    rds = saveRDS(x, file, ...),
    xlsx = openxlsx::write.xlsx(x, file, ...),
    stop("Unknown file type: ", type)
  )
  file
}

write_and_return_file.character <- function(x, file, ...) {
  ensure_directory(file)
  writeLines(x, file, ...)
  file
}

write_and_return_file.ggplot <- function(x, file, ...) {
  ensure_directory(file)
  ggplot2::ggsave(file, plot = x, ...)
  file
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

# align with MAFFT-ginsi
align_mafft_ginsi <- function(seqs, out_file, ncpu, log = "") {
  seqs_file <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(seqs, seqs_file)
  on.exit(unlink(seqs_file))
  args <- c("--globalpair", "--maxiterate", "1000", "--thread", ncpu, seqs_file)
  stopifnot(system2("mafft", args = args, stdout = out_file, stderr = log) == 0)
  out_file
}

iqtree_extensions <- c("treefile", "iqtree", "mldist", "contree", "log", "model.gz", "splits.nex")

# ML tree with IQ-TREE
iqtree <- function(aln, ncpu, log = "") {
  args <- c(
    "-seed", .Random.seed[1],
    "-m", "MFP", # model finder plus
    "-B", "1000" # 1000 ultrafast bootstrap
  )
  # IQ-TREE will auto-resume or just return the old result if called again
  # with the same output prefix.  Generate an output prefix which is unique to
  # the alignment and args (but not necessarily the number of cores)

  alnhash <- tools::md5sum(aln)
  commandhash <- digest::digest(c(alnhash, args))
  alndir <- dirname(aln)
  alnext <- sub(".+\\.", "", basename(aln))
  tempaln <- file.path(alndir, paste(commandhash, alnext, sep = "."))
  file.symlink(basename(aln), tempaln)
  on.exit(unlink(tempaln))

  # now add the alignment name and the cpu specification
  args <- c(
    args,
  "-s", tempaln,
  "-T", ncpu # max number of CPUs, it might not use them efficiently
  )

  # check that all the output files exist; if not force IQ-TREE to run
  tempfiles <- paste(tempaln, iqtree_extensions, sep = ".")
  if (!all(file.exists(tempfiles))) args <- c(args, "-redo")

  stopifnot(
    system2("iqtree", args = args, stdout = log) == 0 ||
      any(grepl("indicates that a previous run successfully finished", readLines(log)))
  )

  # make links from the temp files to the output files.
  outfiles <- paste(aln, iqtree_extensions, sep = ".")
  unlink(outfiles)
  file.link(tempfiles, outfiles)
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

parse_clusters <- function(clusters, singletons, name = "seq_id") {
  tibble::enframe(unname(clusters), name = "cluster", value = "seq_id") %>%
    dplyr::mutate(
      cluster = formatC(cluster, flag = "0",
                        width = ceiling(log10(max(cluster)))),
      seq_id = trimws(seq_id)
    ) %>%
    tidyr::separate_rows(seq_id) %>%
    dplyr::bind_rows(singletons) %>%
    dplyr::rename(!!name := seq_id)
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
    draw_cluster(poly_clusters$mrca, offset = offset, barsize = 2.5, barcolor = scales::alpha("red", 0.5)),
    draw_cluster(mono_clusters$mrca, offset = offset, barsize = 2)
  )
}
cluster_annotation <- function(label, x, y = -0.20, width = 0.05) {
  ggplot2::annotation_custom(
    grid::textGrob(label = label, rot = 90, hjust = 1, vjust = 0.5,
                   gp = grid::gpar(fontsize = 6)),
    xmin = x, xmax = x + width, ymin = y, ymax = y)
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

#### Taxonomy plot ####
# Here's the somewhat modified function I used to make the plots in my previous paper
# Plot distribution between different taxa
taxon_plot <- function(
  .data, # the data
  rank, # what "rank" (column in the data) we are plotting
  ..., # conditions to filter which rows we care about
  y = reads, # should be reads or OTUs
  x = Type, # what groups do we want to divide by
  weight = if ("weight" %in% names(.data)) "weight" else "1", # column to weight the read and ASV counts by
  cutoff = NULL, # groups which represent less than this fraction in all types are grouped together as "other"
  cutoff_type = c("single", "either", "both"),
  data_only = FALSE # just return the data
) {
  rank <- rlang::enquo(rank)
  rank_label <- rlang::as_label(rank)
  y <- rlang::enquo(y)
  x <- rlang::enquo(x)
  weight <- rlang::parse_expr(weight)
  ranks <- c("kingdom", "phylum", "class", "order", "family", "genus")
  cutoff_type <- match.arg(cutoff_type)
  .data <- .data %>%
    dplyr::group_by(!!x) %>%
    dplyr::mutate(
      OTUs = dplyr::n_distinct(OTU),
      reads = reads/sum(reads * !!weight)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(...)
  if (rank_label %in% ranks) {
    .data <- .data %>%
      dplyr::arrange_at(ranks) %>%
      dplyr::mutate_at(
        ranks,
        ~ factor(
          .,
          levels = c(NA, "other",
                     purrr::discard(unique(as.character(.)), is.na)),
          exclude = "NULL"
        )
      )
  } else if (!is.factor(dplyr::pull(.data, !!rank))) {
    .data <- .data %>%
      dplyr::mutate_at(
        rank_label,
        ~ factor(
          .,
          levels = c(NA, "other",
                     purrr::discard(unique(as.character(.)), is.na)),
          exclude = "NULL"
        )
      )
  }
  .data <- .data %>%
    dplyr::group_by(!!x, !!rank) %>%
    dplyr::summarize(reads = sum(!!weight * reads), OTUs = sum(unique(data.frame(OTU, w = !!weight))$w)/max(OTUs)) %>%
    dplyr::ungroup()

  if (data_only) return(.data)

  if (!is.null(cutoff)) {
    prelevels <- levels(dplyr::pull(.data, !!rank))
    cutoff_fun <- switch(
      cutoff_type,
      single = function(x) all(dplyr::pull(x, !!y) < cutoff),
      either = function(x) all(x$reads < cutoff, x$OTUs < cutoff),
      both = function(x) all(x$reads < cutoff) | all(x$OTUs < cutoff)
    )
    .data <- dplyr::group_by(.data, !!rank) %>%
      dplyr::group_map(

        ~ if (cutoff_fun(.x)) {
          dplyr::mutate(
            .x,
            !!rank := factor(
              "other",
              levels = levels(!!rank)
            ),
            exclude = NULL)
        } else {
          .x
        },
        .keep = TRUE
      ) %>%
      dplyr::bind_rows() %>%
      dplyr::group_by(!!x, !!rank) %>%
      dplyr::summarize(reads = sum(reads), OTUs = sum(OTUs)) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(!!rank := factor(!!rank, levels = prelevels, exclude = NULL))
  }

  .data <- dplyr::mutate(.data, !!rank := forcats::fct_drop(!!rank))
  vals <- levels(dplyr::pull(.data, !!rank))
  # if ("other" %in% vals) vals <- c("other", vals) %>% magrittr::extract(!duplicated(.))
  # if (any(is.na(vals))) vals <- c(NA, vals) %>% magrittr::extract(!duplicated(.))

  if (rank_label == stringr::str_to_lower(rank_label)) rank_label <- stringr::str_to_title(rank_label)
  y_label <- rlang::as_label(y)
  if (y_label == stringr::str_to_lower(y_label)) y_label <- stringr::str_to_title(y_label)
  ggplot(.data, aes(x = !!x, y = !!y, fill = !!rank)) +
    geom_bar(position = position_stack(reverse = TRUE),
             stat = "identity", color = "white", size = 0.2) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          strip.background = element_blank(),
          panel.spacing = unit(3, "pt")) +
    scale_fill_discrete(
      breaks = vals,
      labels = tidyr::replace_na(as.character(vals), "unidentified"),
      name = rank_label,
      guide = guide_legend(ncol = 4, byrow = TRUE)
    ) +
    ylab(paste("Fraction of", y_label)) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "bottom") +
    xlab(NULL)

}

`%||%` <- function(x, ifna) {
  ifelse(is.na(x), ifna, x)
}

tree_depth <- function(phylo) {
  n <- length(phylo$edge.length)
  depth <- numeric(n)
  for (i in seq_len(n)) {
    depth[i] <- phylo$edge.length[i] + depth[match(phylo$edge[i, 1], phylo$edge[,2])] %||% 0
  }
  max(depth)
}
