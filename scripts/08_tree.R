# Compare results from different denoising pipelines by making a tree
# For the tree, use only 5.8S and the full LSU regions
# these can be aligned separately and concatenated.

library(targets)
library(tarchetypes)
library(rlang)

#### Combined table ####
align_targets <- tar_map(
  values = list(
    hash = c("5_8S_hash", "LSU_hash", "ITS2_hash", "ITS_hash"),
    region = c("5_8S", "LSU", "ITS2", "ITS")
  ),
  # get all the unique 5.8S and LSU sequences
  tar_combine(
    allseqs,
    positions_targets$regions,
    command =
      purrr::map_dfr(list(!!!.x), dplyr::select, hash, region) %>%
      unique() %>%
      dplyr::filter(complete.cases(.)) %>%
      tibble::deframe() %>%
      chartr(old = "T", new = "U") %>%
      Biostrings::RNAStringSet()
  ),
  #### Alignment ####
  # align them independently
  # these are pretty quick alignments (~1 hour)
  # they are a bit faster because we are only aligning the unique sequences
  tar_target(
    name = align,
    DECIPHER::AlignSeqs(allseqs, processors = 8) %>%
      # add an extra sequence to each one for "missing"
      # it is just gaps, with the same length as the alignment
      magrittr::inset(
        i = "missing",
        value = Biostrings::RNAStringSet(strrep("-", Biostrings::width(.[1])))
      )
  ),
  names = region
)

# we don't need to align ITS
align_targets$align <- purrr::discard(
  align_targets$align,
  ~ grepl("ITS", .$settings$name)
)


concat_targets <- tar_plan(

  tar_combine(
    hash_key,
    positions_targets$hash_key,
    command = dplyr::bind_rows(!!!.x) %>% unique()
  ),

  #### Concatenation ####
  # now find all unique pairs of 5.8S and LSU
  tar_combine(
    all_both,
    positions_targets$regions,
    command =
      purrr::map_dfr(list(!!!.x), dplyr::select, "5_8S_hash", "LSU_hash") %>%
      unique()
  ),
  # paste together the 5.8S and LSU sequences for each.
  align_both = paste(
    align_5_8S[all_both$`5_8S_hash`],
    align_LSU[all_both$LSU_hash],
    sep = ""
  ) %>%
    set_names(paste(all_both$`5_8S_hash`, all_both$LSU_hash, sep = "_")) %>%
    # make it DNA instead of RNA for fastree
    chartr(old = "U", new = "T") %>%
    Biostrings::DNAStringSet(),

  # remove columns with at least 90% gaps
  align_degap = Biostrings::DNAMultipleAlignment(align_both) %>%
    Biostrings::maskGaps(min.fraction = 0.9, min.block.width = 1) %>%
    as("DNAStringSet"),

  tar_map(
    values = list(
      aln = rlang::syms(c("align_both", "align_degap")),
      file = c("comparealn", "comparealn.degap"),
      id = c("both", "degap")
    ),
    tar_file(
      write_aln,
      {
        if (!dir.exists(comparedir)) dir.create(comparedir)
        alnfile <- file.path(comparedir, sprintf("%s.fasta", file))
        Biostrings::writeXStringSet(aln, alnfile)
        alnfile
      }
    ),
    names = id
  ),

  #### ML Tree ####
  # make a tree with fasttree
  # takes about 5 min
  tar_file(
    tree_file,
    fasttree(
      aln_file = write_aln_degap,
      out_file = file.path(comparedir, "comparealn.degap.tree"),
      constraints = tax_constraints
    )
  ),

  tree_raw = treeio::read.newick(tree_file),
  tree_rooted = root_with_kingdoms(tree_raw, kingdoms, c("Euglenozoa", "Heterolobosa")),
  tree_fungi = dplyr::filter(kingdoms, taxon == "Fungi")$label %>%
    ape::getMRCA(tree_rooted, .) %>%
    ape::extract.clade(tree_rooted, .),
  tree_animals = dplyr::filter(kingdoms, taxon == "Metazoa")$label %>%
    ape::getMRCA(tree_rooted, .) %>%
    ape::extract.clade(tree_rooted, .),
  tree_plants = dplyr::filter(phyla, taxon == "Streptophyta")$label %>%
    ape::getMRCA(tree_rooted, .) %>%
    ape::extract.clade(tree_rooted, .),
  tree_protists = tree_rooted %>%
    ape::drop.tip(tree_fungi$tip.label) %>%
    ape::drop.tip(tree_animals$tip.label) %>%
    ape::drop.tip(tree_plants$tip.label)
)

#### Community matrix ####
# make a "community matrix" for the different pipelines
reads_targets <- tar_map(
  tidyr::crossing(
    m = tibble::tibble(
      table = rlang::syms(c("ampliseq_table", "vs_table", "sl_table")),
      regions = rlang::syms(c("regions_as", "regions_vs", "regions_sl")),
      id = c("ampliseq", "vsearch", "single_link")
    ),
    r = tibble::tibble(
      hashcols = list(c("5_8S_hash", "LSU_hash"), "ITS2_hash", "ITS_hash"),
      id2 = c("concat", "ITS2", "ITS")
    )
  ) %>% tidyr::unpack(c("m", "r")),
  tar_fst_tbl(
    reads,
    tibble::column_to_rownames(table, "OTU") %>%
      rowSums() %>%
      tibble::enframe(name = "seq_id", value = "n") %>%
      dplyr::left_join(regions, by = "seq_id") %>%
      dplyr::mutate_at(hashcols, tidyr::replace_na, "missing") %>%
      dplyr::transmute(
        seq_id = glue::glue(paste0("{`", hashcols, "`}", collapse = "_")),
        n = n/sum(n)
      ) %>%
      dplyr::rename(!!id := n) %>%
      dplyr::group_by(seq_id) %>%
      dplyr::summarize_all(sum),
    tidy_eval = FALSE
  ),
  names = c(id, id2)
)

phyloseq_targets <- tar_plan(

  tar_combine(
    otu_tab,
    purrr::keep(reads_targets$reads, ~ endsWith(.$settings$name, "concat")),
    command =
      purrr::reduce(list(!!!.x), dplyr::full_join, by = "seq_id") %>%
      dplyr::mutate_if(is.numeric, tidyr::replace_na, 0) %>%
      tibble::column_to_rownames("seq_id") %>%
      phyloseq::otu_table(taxa_are_rows = TRUE)
  ),

  #### phyloseq object and UniFrac distances ####
  tax_table =
    dplyr::select(tax_all, full_hash, rank, taxon) %>%
    dplyr::left_join(
      dplyr::transmute(
        hash_key,
        full_hash = full_hash,
        label = stringr::str_c(`5_8S_hash`, LSU_hash, sep = "_")
      ),
      by = "full_hash"
    ) %>%
    phylotax::make_taxon_labels(abbrev = taxon_abbrevs) %>%
    dplyr::rename(taxon_label = new) %>%
    dplyr::left_join(tibble::tibble(old = tree_rooted$tip.label), by = "old") %>%
    tidyr::replace_na(list(taxon_label = "")) %>%
    tibble::column_to_rownames("old") %>%
    as.matrix() %>%
    phyloseq::tax_table(),
  tar_map(
    values = list(
      tree = rlang::syms(c("tree_rooted", "tree_fungi")),
      id = c("alleuks", "fungi")
    ),
    tar_target(
      physeq,
      phyloseq::phyloseq(tree, otu_tab, tax_table)
    ),
    names = id
  ),
  # physeq_glom = phyloseq::tip_glom(physeq, h = 0.01),

  tar_map(
    values = tidyr::crossing(
      weight = c(TRUE, FALSE),
      tidyr::nesting(
        physeq = rlang::syms(c("physeq_alleuks", "physeq_fungi")),
        id = c("alleuks", "fungi")
      )
    ) %>%
      dplyr::mutate(id = paste0(id, "_", ifelse(weight, "", "un"), "weighted")),
    tar_target(
      unifrac,
      phyloseq::UniFrac(physeq, weighted = weight)
    ),
    names = id
  ),

  #### Figure ####

  tar_target(
    fungi_cluster_geoms,
    c(
      draw_clusters(
        clusters = its2_cluster_90,
        singletons = its2_precluster_singletons,
        hash_key = hash_key,
        physeq = physeq_fungi,
        offset = 0.15
      ),
      draw_clusters(
        clusters = its2_cluster_97,
        singletons = its2_precluster_singletons,
        hash_key = hash_key,
        physeq = physeq_fungi,
        offset = 0.20
      ),
      draw_clusters(
        clusters = its2_cluster_99,
        singletons = its2_precluster_singletons,
        hash_key = hash_key,
        physeq = physeq_fungi,
        offset = 0.25
      )
    ),
    packages = c("rlang", "ggtree")
  ),

  tar_map(
    values = list(
      physeq_obj = rlang::syms(c("physeq_alleuks", "physeq_fungi")),
      tree_height = c(300, 75),
      tip_offset = c(0.4, 0.3),
      id = c("alleuks", "fungi")
    ),
    tar_qs(
      tree_fig,
      (ggtree(physeq_obj) +
         theme_tree2() +
         geom_tiplab(label = "", align = TRUE) +
         geom_tiplab(
           aes(label = taxon_label),
           align = TRUE, offset = tip_offset, linetype = NULL, size = 2.5, family = "mono") +
         xlim(0, 6) +
         scale_y_continuous(expand = expansion(add = c(10, 0.6)))) %>%
        gheatmap(
          log10(as(phyloseq::otu_table(physeq_obj), "matrix")),
          colnames_angle = 90,
          hjust = 1,
          width = 0.1,
          legend_title = "log10(read abundance)",
          font.size = 2
        ) +
        ggplot2::theme(legend.position = "bottom"),
      packages = c("ggplot2", "ggtree", "rlang"),
      error = "workspace"
    ),
    tar_file(
      tree_fig_file,
      sprintf("%s/treemap_%s.pdf", figdir, id) %T>%
        ggplot2::ggsave(., plot = tree_fig, device = "pdf", width = 12,
                        height = tree_height, limitsize = FALSE)
    ),
    names = id
  ),
  tar_file(
    tree_cluster_fig_file,
    sprintf("%s/tree_clusters_fungi.pdf", figdir) %T>%
      ggplot2::ggsave(
        .,
        plot = tree_fig_fungi +
          fungi_cluster_geoms +
          cluster_annotation(label = "SH 0.90", x = 1.425) +
          cluster_annotation(label = "SH 0.97", x = 1.475) +
          cluster_annotation(label = "SH 0.99", x = 1.525),
        device = "pdf", width = 12, height = 75, limitsize = FALSE)
  )
)

tree_targets <- c(
  align_targets,
  concat_targets,
  reads_targets,
  phyloseq_targets
)

# colSums(phyloseq::otu_table(physeq_glom) > 0)
# colSums(phyloseq::otu_table(physeq_glom) > 1/30000) # no global singletons


# How many unique of each?
# dplyr::n_distinct(ampliseq_regions$`5_8S`)
# dplyr::n_distinct(ampliseq_regions$LSU)

#### Tzara Clusters ####
# now load the tzara data
# tzara_sets <- c("quiver_nosingle", "quiver_single", "arrow_nosingle", "arrow_single")
# tzara_regions <- list()
# for (s in tzara_sets) {
#     tzara_regions[[s]] <- list()
#     for (r in regions) {
#         tzara_regions[[s]][[r]] <- Biostrings::readDNAStringSet(
#             paste0("processReads/tzara/", s, "/ASV_", r, ".fasta")
#         )
#         tzara_regions[[s]][[r]] <- tibble::enframe(
#             as.character(tzara_regions[[s]][[r]]),
#             name = "seq_id",
#             value = r
#         )
#     }
#     tzara_regions[[s]] <-
#         purrr::reduce(tzara_regions[[s]], dplyr::full_join, by = "seq_id")
#     tzara_regions[[s]] <- dplyr::transmute(
#         tzara_regions[[s]],
#         seq_id = seq_id,
#         `5_8S` = `5_8S`,
#         `5_8S_hash` = tzara::seqhash(`5_8S`),
#         LSU = stringr::str_c(LSU1, V2, LSU2, V3, LSU3, V4, LSU4),
#         LSU_hash = tzara::seqhash(LSU)
#     )
#     tzara_regions[[s]] <- dplyr::mutate_at(
#         tzara_regions[[s]],
#         dplyr::vars(dplyr::ends_with("_hash")),
#         tidyr::replace_na,
#         replace = "missing"
#     )
# }

# How many unique of each?
# dplyr::n_distinct(sl_regions$`5_8S`)
# dplyr::n_distinct(sl_regions$LSU)

#### VSEARCH clusters ####

# How many unique of each?
# dplyr::n_distinct(vs_regions$`5_8S`)
# dplyr::n_distinct(vs_regions$LSU)

#### LAA ####

# laa_seqs <- Biostrings::readDNAStringSet(
#     here::here("process", "pb_363_subreads.demux.sieve.swarm.laa.select.fastq.gz"),
#     format = "fastq"
# )
#
# laa_positions <- LSUx::lsux(
#     seq = laa_seqs,
#     cm_32S = cm_32S_trunc,
#     ITS1 = TRUE,
#     cpu = 8,
#     # allow 2 Gb ram (per process)
#     mxsize = 2048
# )
#
# # Just cut out 5.8S and LSU
# laa_regions <- purrr::map2(
#     .x = c("5_8S", "LSU1"),
#     .y = c("5_8S", "LSU4"),
#     tzara::extract_region,
#     seq = as.character(laa_seqs),
#     positions = laa_positions
# )
# laa_regions <- purrr::map2(laa_regions,
#                           c("5_8S", "LSU"),
#                           tibble::enframe,
#                           name = "seq_id")
# laa_regions <- purrr::reduce(
#     laa_regions,
#     dplyr::full_join,
#     by = "seq_id"
# )
# laa_regions <- dplyr::mutate(
#     laa_regions,
#     `5_8S_hash` = tzara::seqhash(`5_8S`),
#     LSU_hash = tzara::seqhash(LSU)
# )
#
# laa_regions <- dplyr::mutate_at(
#     laa_regions,
#     dplyr::vars(dplyr::ends_with("_hash")),
#     tidyr::replace_na,
#     replace = "missing"
# )
#
# # How many unique of each?
# dplyr::n_distinct(laa_regions$`5_8S`)
# dplyr::n_distinct(laa_regions$LSU)
