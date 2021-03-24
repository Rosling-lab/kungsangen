# Compare results from different denoising pipelines by making a tree
# For the tree, use only 5.8S and the full LSU regions
# these can be aligned separately and concatenated.

library(targets)
library(tarchetypes)

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
                id = c("as", "sl", "vs")),
  tar_target(
    positions,
    LSUx::lsux(
      seq = seq,
      cm_32S = cm_32S_trunc,
      ITS1 = TRUE,
      cpu = 8,
      # allow 2 Gb ram (per process)
      mxsize = 2048
    )
  ),

  # Just cut out 5.8S and LSU
  tar_target(
    regions,
    purrr::map2(
      .x = c("5_8S", "LSU1"),
      .y = c("5_8S", "LSU4"),
      tzara::extract_region,
      seq = seq,
      positions = positions
    ) %>%
      purrr::map2(c("5_8S", "LSU"), tibble::enframe, name = "seq_id") %>%
      purrr::reduce(dplyr::full_join, by = "seq_id") %>%
      dplyr::mutate(
        `5_8S_hash` = tzara::seqhash(`5_8S`),
        LSU_hash = tzara::seqhash(LSU)
      ) %>%
      dplyr::mutate_at(
        dplyr::vars(dplyr::ends_with("_hash")),
        tidyr::replace_na,
        replace = "missing"
      )
  ),
  names = id
)


#### Combined table ####
align_targets <- tar_map(
  values = list(hash = c("5_8S_hash", "LSU_hash"), region = c("5_8S", "LSU")),
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

concat_targets <- tar_plan(
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

  # output the alignment
  tar_file(
    write_comparealn,
    {
      comparedir <- file.path("processReads", "compare")
      if (!dir.exists(comparedir)) dir.create(comparedir)
      alnfile <- file.path(comparedir, "comparealn.fasta.gz")
      Biostrings::writeXStringSet(align_both, alnfile, compress = TRUE)
      alnfile
    }
  ),

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
        comparedir <- file.path("processReads", "compare")
        if (!dir.exists(comparedir)) dir.create(comparedir)
        alnfile <- file.path(comparedir, sprintf("%s.fasta.gz", file))
        Biostrings::writeXStringSet(aln, alnfile, compress = TRUE)
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
    file.path(comparedir, "comparealn.degap.tree") %T>%
    system2(
      command = "fasttree",
      args = c("-nt", "-gtr", write_aln_degap),
      stdout = .
    )
  ),

  tree = treeio::read.newick(tree_file) %>%
    phangorn::midpoint()
)

#### Community matrix ####
# make a "community matrix" for the different pipelines
reads_targets <- tar_map(
  list(
    table = rlang::syms(c("ampliseq_table", "vs_table", "sl_table")),
    regions = rlang::syms(c("regions_as", "regions_vs", "regions_sl")),
    id = c("ampliseq", "vsearch", "single_link")
  ),
  tar_fst_tbl(
    reads,
    tibble::column_to_rownames(table, "OTU") %>%
      rowSums() %>%
      tibble::enframe(name = "seq_id", value = "n") %>%
      dplyr::left_join(regions, by = "seq_id") %>%
      dplyr::mutate_at(c("5_8S_hash", "LSU_hash"), tidyr::replace_na, "missing") %>%
      dplyr::transmute(
        seq_id = paste(`5_8S_hash`, LSU_hash, sep = "_"),
        n = n/sum(n)
      ) %>%
      dplyr::rename(!!id := n) %>%
      dplyr::group_by(seq_id) %>%
      dplyr::summarize_all(sum),
    tidy_eval = FALSE
  ),
  names = id
)

phyloseq_targets <- tar_plan(

  tar_combine(
    otu_tab,
    reads_targets$reads,
    command =
      purrr::reduce(list(!!!.x), dplyr::full_join, by = "seq_id") %>%
      dplyr::mutate_if(is.numeric, tidyr::replace_na, 0) %>%
      tibble::column_to_rownames("seq_id") %>%
      phyloseq::otu_table(taxa_are_rows = TRUE)
  ),

  #### phyloseq object and UniFrac distances ####
  physeq = phyloseq::phyloseq(tree, otu_tab),
  physeq_glom = phyloseq::tip_glom(physeq, h = 0.01),

  tar_map(
    values = list(weight = c(TRUE, FALSE), id = c("unweighted", "weighted")),
    tar_target(
      unifrac,
      phyloseq::UniFrac(physeq, weighted = weight)
    ),
    names = id
  ),

  #### Figure ####
  tar_map(
    values = list(
      physeq_obj = rlang::syms(c("physeq", "physeq_glom")),
      id = c("raw", "glom")
    ),
    tar_qs(
      tree_fig,
      (ggtree::ggtree(physeq_obj) +
         ggtree::theme_tree2() +
         ggtree::geom_tiplab(label = "", align = TRUE)) %>%
        ggtree::gheatmap(log(as(phyloseq::otu_table(physeq_obj), "matrix")),
                         colnames_angle = 90,
                         hjust = 1,
                         width = 0.2,
                         legend_title = "log10(read abundance)")
    ),
    tar_file(
      tree_fig_file,
      sprintf("processReads/compare/treemap_%s.pdf", id) %T>%
        ggplot2::ggsave(., plot = tree_fig, device = "pdf", width = 8,
                        height = 300, limitsize = FALSE)
    ),
    names = id
  )
)

compare_tree_targets <- list(
  pre_positions_targets,
  positions_targets,
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
