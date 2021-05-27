# Compare results from different denoising pipelines by making a tree
# For the tree, use only 5.8S and the full LSU regions
# these can be aligned separately and concatenated.

library(targets)
library(tarchetypes)
library(rlang)

#### Combined table ####
align_plan <- tar_map(
  values = list(
    hash = c("5_8S_hash", "LSU_hash", "ITS2_hash", "ITS_hash"),
    region = c("5_8S", "LSU", "ITS2", "ITS")
  ),
  #### allseqs_{region} ####
  # get all the unique 5.8S and LSU sequences
  tar_combine(
    allseqs,
    positions_plan$regions,
    command =
      purrr::map_dfr(list(!!!.x), dplyr::select, hash, region) %>%
      unique() %>%
      dplyr::filter(complete.cases(.)) %>%
      tibble::deframe() %>%
      chartr(old = "T", new = "U") %>%
      Biostrings::RNAStringSet()
  ),
  #### *Alignment* ####
  # align them independently
  # these are pretty quick alignments (~1 hour)
  # they are a bit faster because we are only aligning the unique sequences
  #### align_{region} ####
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
align_plan$align <- purrr::discard(
  align_plan$align,
  ~ grepl("ITS", .$settings$name)
)


concat_plan <- tar_plan(
  #### hash_key ####
  tar_combine(
    hash_key,
    positions_plan$hash_key,
    command = dplyr::bind_rows(!!!.x) %>% unique()
  ),

  #### *Concatenation* ####
  #### align_key ####
  # now find all unique pairs of 5.8S and LSU
  align_key = dplyr::select(hash_key, "label", "5_8S_hash", "LSU_hash") %>%
    unique(),
  #### align_both ####
  # paste together the 5.8S and LSU sequences for each.
  align_both = paste(
    align_5_8S[align_key$`5_8S_hash`],
    trim_LSU_intron(align_LSU)[align_key$LSU_hash],
    sep = ""
  ) %>%
    set_names(align_key$label) %>%
    # make it DNA instead of RNA for fastree
    chartr(old = "U", new = "T") %>%
    Biostrings::DNAStringSet(),

  #### align_degap ####
  # remove columns with at least 90% gaps
  align_degap = Biostrings::DNAMultipleAlignment(align_both) %>%
    Biostrings::maskGaps(min.fraction = 0.9, min.block.width = 1) %>%
    as("DNAStringSet"),

  tar_map(
    values = list(
      aln = rlang::syms(c("align_both", "align_degap")),
      file = c("comparealn", "comparealn.degap"),
      degap_id = c("both", "degap")
    ),
    #### write_aln_{degap_id} ####
    tar_file(
      write_aln,
      write_and_return_file(aln, file.path(comparedir, sprintf("%s.fasta", file)))
    ),
    names = degap_id
  ),

  #### *ML Tree* ####
  #### tree_file ####
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

  #### *separate out Fungi and protist trees* ####
  #### tree_raw ####
  tree_raw = treeio::read.newick(tree_file),
  #### tree_rooted ####
  tree_rooted = root_with_kingdoms(tree_raw, kingdoms, c("Euglenozoa", "Heterolobosa")),
  #### tree_fungi ####
  tree_fungi = dplyr::filter(kingdoms, taxon == "Fungi")$label %>%
    ape::getMRCA(tree_rooted, .) %>%
    ape::extract.clade(tree_rooted, .),
  #### tree_animals ####
  tree_animals = dplyr::filter(kingdoms, taxon == "Metazoa")$label %>%
    ape::getMRCA(tree_rooted, .) %>%
    ape::extract.clade(tree_rooted, .),
  #### tree_plants ####
  tree_plants = dplyr::filter(phyla, taxon == "Streptophyta")$label %>%
    ape::getMRCA(tree_rooted, .) %>%
    ape::extract.clade(tree_rooted, .),
  #### tree_protists ####
  tree_protists = tree_rooted %>%
    ape::drop.tip(tree_fungi$tip.label) %>%
    ape::drop.tip(tree_animals$tip.label) %>%
    ape::drop.tip(tree_plants$tip.label),

  #### *Realign and retree for Fungi* ####
  # pick the most abundant nonfungal obazoan to use as an outgroup for Fungi
  #### fungi_outgroup ####
  fungi_outgroup = phyloseq::subset_taxa(
    physeq_alleuks,
    grepl("Ichthyosporia|Choanoflagello|Meta-", taxon_label)
  ) %>%
    phyloseq::otu_table() %>%
    rowSums() %>%
    which.max() %>%
    names() %>%
    strsplit("_") %>%
    unlist() %>%
    set_names(c("5_8S", "LSU")),

  #### regions_fungi ####
  # pick out the fungi
  regions_fungi = tibble::tibble(label = tree_fungi$tip.label) %>%
    tidyr::separate(label, c("5_8S", "LSU")) %>%
    as.list(),

  tar_map(
    values = tibble::tibble(
      region = c("5_8S", "LSU"),
      allseqs = paste0("allseqs_", region) %>% rlang::syms()
    ),
    names = region,
    #### realign_{region} ####
    tar_file(
      realign,
      align_mafft_ginsi(
        seqs = allseqs[c(fungi_outgroup[[region]], unique(regions_fungi[[region]]))],
        out_file = file.path(comparedir, sprintf("fungi_realign_%s.fasta", region)),
        ncpu = local_cpus(),
        log = file.path("logs", sprintf("fungi_realign_%s.log", region))
      )
    ),
    #### realign_copies_{region} ####
    tar_target(
      realign_copies,
      Biostrings::readRNAStringSet(realign)[
        c(fungi_outgroup[[region]], regions_fungi[[region]])
        ]
    ),
    unlist = TRUE
  ),
  #### reconcat ####
  tar_file(
    reconcat,
    paste(
      realign_copies_5_8S,
      trim_LSU_intron(realign_copies_LSU),
      sep = ""
    ) %>%
    set_names(paste(names(realign_copies_5_8S), names(realign_copies_LSU), sep = "_")) %>%
    chartr(old = "Uu", new = "Tt") %>%
    Biostrings::DNAStringSet() %>%
      write_and_return_file(file.path(comparedir, "fungi_realign.fasta"))
  ),
  #### iqtree_fungi ####
  tar_file(
    iqtree_fungi,
    iqtree(reconcat, local_cpus())
  ),
  #### tree_fungi_new ####
  tar_target(
    tree_fungi_new,
    ape::read.tree(iqtree_fungi[1]) %>%
      ape::drop.tip(paste(fungi_outgroup, collapse = "_"))
  )
)

#### *Community matrix* ####
# make a "community matrix" for the different pipelines
reads_plan <- tar_map(
  tidyr::crossing(
    m = tibble::tibble(
      table = rlang::syms(c("ampliseq_table", "table_vs", "table_sl")),
      regions = rlang::syms(c("regions_as", "regions_vs", "regions_sl")),
      otuid = c("ampliseq", "vsearch", "single_link")
    ),
    r = tibble::tibble(
      hashcols = list(c("5_8S_hash", "LSU_hash"), "ITS2_hash", "ITS_hash", "full_hash"),
      regionid = c("concat", "ITS2", "ITS", "full")
    )
  ) %>% tidyr::unpack(c("m", "r")),
  #### reads_{otuid}_{regionid} ####
  tar_fst_tbl(
    reads,
    tibble::column_to_rownames(table, "OTU") %>%
      rowSums() %>%
      tibble::enframe(name = "seq_id", value = "n") %>%
      dplyr::left_join(regions, by = "seq_id") %>%
      dplyr::mutate_at(hashcols, tidyr::replace_na, "missing") %>%
      dplyr::transmute(
        seq_id = glue::glue(paste0("{`", hashcols, "`}", collapse = "_")),
        n = n
      ) %>%
      dplyr::rename(!!otuid := n) %>%
      dplyr::group_by(seq_id) %>%
      dplyr::summarize_all(sum),
    tidy_eval = FALSE
  ),
  #### abundance_{otuid}_{regionid} ####
  tar_fst_tbl(
    abundance,
    dplyr::mutate_if(reads, is.numeric, ~./sum(.))
  ),
  names = c(otuid, regionid)
)

phyloseq_plan <- tar_plan(
  #### reads_tab ####
  tar_combine(
    reads_tab,
    purrr::keep(reads_plan$reads, ~ endsWith(.$settings$name, "ITS2")),
    command =
      purrr::reduce(list(!!!.x), dplyr::full_join, by = "seq_id") %>%
      dplyr::mutate_if(is.numeric, tidyr::replace_na, 0L)
  ),

  #### otu_tab ####
  tar_combine(
    otu_tab,
    purrr::keep(reads_plan$abundance, ~ endsWith(.$settings$name, "concat")),
    command =
      purrr::reduce(list(!!!.x), dplyr::full_join, by = "seq_id") %>%
      dplyr::mutate_if(is.numeric, tidyr::replace_na, 0) %>%
      tibble::column_to_rownames("seq_id") %>%
      phyloseq::otu_table(taxa_are_rows = TRUE)
  ),

  #### *phyloseq object and UniFrac distances* ####
  #### tax_table ####
  tax_table =
    dplyr::select(tax_all, full_hash, rank, taxon) %>%
    dplyr::left_join(
      dplyr::select(hash_key, full_hash, label),
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
      tree = rlang::syms(c("tree_rooted", "tree_fungi_new")),
      group = c("alleuks", "fungi")
    ),
    #### physeq_{group} ####
    tar_target(
      physeq,
      phyloseq::phyloseq(tree, otu_tab, tax_table)
    ),
    names = group
  ),
  # physeq_glom = phyloseq::tip_glom(physeq, h = 0.01),

  tar_map(
    values = tidyr::crossing(
      weight = c(TRUE, FALSE),
      tidyr::nesting(
        physeq = rlang::syms(c("physeq_alleuks", "physeq_fungi")),
        group = c("alleuks", "fungi")
      )
    ) %>%
      dplyr::mutate(group = paste0(group, "_", ifelse(weight, "", "un"), "weighted")),
    #### unifrac_{group}_{weighted} ####
    tar_target(
      unifrac,
      phyloseq::UniFrac(physeq, weighted = weight)
    ),
    names = group
  ),

  #### Figure ####
  tar_map(
    # make separate figures for fungi and all eukaryotes
    values = list(
      physeq_obj = rlang::syms(c("physeq_alleuks", "physeq_fungi")),
      tree_height = c(150, 75),
      tip_offset = c(0.4, 0.3),
      group = c("alleuks", "fungi")
    ),
    #### tree_fig_{group} ####
    # just the tree without the clusters
    tar_qs(
      tree_fig,
      (
        ggtree(physeq_obj) +
          theme_tree2() +
          geom_tiplab(label = "", align = TRUE) +
          scale_y_continuous(expand = expansion(add = c(5, 0.6)))

      ) %>%
        {
          . +
            geom_tiplab(
              aes(label = taxon_label),
              align = TRUE,
              offset = max(.$data$x) * 0.19,
              linetype = NULL,
              size = 2.5,
              family = "mono"
            ) +
            scale_x_continuous(
              expand = expansion(mult = c(0.01, 1.5)),
              guide = NULL
              # breaks = max(.$data$x) * c(1.12, 1.14, 1.16),
              # labels = c("GH 0.90", "SH 0.97", "SH 0.99")
            )
        } %>%
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
    #### tree_fig_file_{group} ####
    tar_file(
      tree_fig_file,
      write_and_return_file(
        tree_fig,
        sprintf("%s/treemap_%s.pdf", figdir, group),
        device = "pdf", width = 12,
        height = tree_height, limitsize = FALSE
      )
    ),
    #### cluster_geoms_{group} ####
    # clusters
    tar_target(
      cluster_geoms,
      c(
        draw_clusters(
          clusters = its2_cluster_90,
          singletons = its2_precluster_singletons,
          hash_key = hash_key,
          physeq = physeq_obj,
          offset = max(tree_fig$data$x) * 0.11
        ),
        draw_clusters(
          clusters = its2_cluster_97,
          singletons = its2_precluster_singletons,
          hash_key = hash_key,
          physeq = physeq_obj,
          offset = max(tree_fig$data$x) * 0.135
        ),
        draw_clusters(
          clusters = its2_cluster_99,
          singletons = its2_precluster_singletons,
          hash_key = hash_key,
          physeq = physeq_obj,
          offset = max(tree_fig$data$x) * 0.16
        )
      ),
      packages = c("rlang", "ggtree")
    ),
    #### tree_cluster_fig_file_{group} ####
    tar_file(
      tree_cluster_fig_file,
      write_and_return_file(
        tree_fig +
          cluster_geoms +
          ggplot2::annotate("text", x = max(tree_fig$data$x) * 1.13, y = 0,
                            label = "GH 90", angle = 90, hjust = 1, size = 2) +
          ggplot2::annotate("text", x = max(tree_fig$data$x) * 1.155, y = 0,
                            label = "SH 97", angle = 90, hjust = 1, size = 2) +
          ggplot2::annotate("text", x = max(tree_fig$data$x) * 1.18, y = 0,
                            label = "SH 99", angle = 90, hjust = 1, size = 2),
        # cluster_annotation(label = "SH 0.90", x = max(tree_fig$data$x) * 1.12) +
        # cluster_annotation(label = "SH 0.97", x = max(tree_fig$data$x) * 1.14) +
        # cluster_annotation(label = "SH 0.99", x = max(tree_fig$data$x) * 1.16),
        sprintf("%s/tree_clusters_%s.pdf", figdir, group),
        device = "pdf", width = 12, height = tree_height, limitsize = FALSE
      )
    ),
    names = group
  )
)

constraint_plan <- list(
  #### constraint_tree_file ####
  constraint_tree_file = tar_file(
    constraint_tree_file,
    file.path("reference", "constraints.tree")
  ),
  #### constraint_tree_plot ####
  constraint_tree_plot = tar_target(
    constraint_tree_plot,
    ape::read.tree(constraint_tree_file) %>%
      purrr::modify_at("node.label", chartr, old = "_", new = " ") %>%
      ggtree(layout = "rectangular", ladderize = FALSE) +
      geom_tiplab(hjust = 1, nudge_y = 0.3, size = 1.7) +
      geom_nodelab(hjust = 1, nudge_y = 0.3, nudge_x = -0.1, size = 1.7),
    packages = "ggtree"
  ),
  constraint_tree_write = tar_map(
    values = plot_type_meta,
    names = ext,
    #### constraint_tree_{ext} ####
    tar_file(
      constraint_tree,
      write_and_return_file(
        constraint_tree_plot,
        file.path(figdir, paste0("constraint_tree.", ext)),
        device = fun, width = 3, height = 2, dpi = 150
      )
    )
  )
)

tree_plan <- c(
  align_plan,
  concat_plan,
  reads_plan,
  phyloseq_plan,
  constraint_plan
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
# dplyr::n_distinct(regions_sl$`5_8S`)
# dplyr::n_distinct(regions_sl$LSU)

#### VSEARCH clusters ####

# How many unique of each?
# dplyr::n_distinct(regions_vs$`5_8S`)
# dplyr::n_distinct(regions_vs$LSU)

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
