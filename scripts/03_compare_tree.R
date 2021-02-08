# Compare results from different denoising pipelines by making a tree
# For the tree, use only 5.8S and the full LSU regions
# these can be aligned separately and concatenated.

#### Ampliseq clusters ####
# First load the ampliseq results and use LSUx to cut out subregions
ampliseq <- readr::read_tsv(
    "processReads/ampliseq/qiime2_ASV_table.tsv",
    col_types = readr::cols(
        .default = readr::col_double(),
        Feature.ID = readr::col_character(),
        Taxon = readr::col_character(),
        sequence = readr::col_character()
    )
)

ampliseq <- tibble::deframe(ampliseq[,c("Feature.ID", "sequence")])

# get the path for the CM which is truncated at the LR5 primer site
# (included in LSUx)
cm_32S_trunc <- system.file(
    file.path("extdata", "fungi_32S_LR5.cm"),
    package = "LSUx"
)

positions <- LSUx::lsux(
    seq = ampliseq,
    cm_32S = cm_32S_trunc,
    ITS1 = TRUE,
    cpu = 8,
    # allow 2 Gb ram (per process)
    mxsize = 2048
)

regions <- unique(positions$region)
regions

# Just cut out 5.8S and LSU
ampliseq_regions <- purrr::map2(
    .x = c("5_8S", "LSU1"),
    .y = c("5_8S", "LSU4"),
    tzara::extract_region,
    seq = ampliseq,
    positions = positions
)
ampliseq_regions <- purrr::map2(ampliseq_regions,
                                c("5_8S", "LSU"),
                                tibble::enframe,
                                name = "seq_id")
ampliseq_regions <- purrr::reduce(
    ampliseq_regions,
    dplyr::full_join,
    by = "seq_id"
)
ampliseq_regions <- dplyr::mutate(
    ampliseq_regions,
    `5_8S_hash` = tzara::seqhash(`5_8S`),
    LSU_hash = tzara::seqhash(LSU)
)
ampliseq_regions <- dplyr::mutate_at(
    ampliseq_regions,
    dplyr::vars(dplyr::ends_with("_hash")),
    tidyr::replace_na,
    replace = "missing"
)

# How many unique of each?
dplyr::n_distinct(ampliseq_regions$`5_8S`)
dplyr::n_distinct(ampliseq_regions$LSU)

#### Tzara Clusters ####
# now load the tzara data
tzara_sets <- c("quiver_nosingle", "quiver_single", "arrow_nosingle", "arrow_single")
tzara_regions <- list()
for (s in tzara_sets) {
    tzara_regions[[s]] <- list()
    for (r in regions) {
        tzara_regions[[s]][[r]] <- Biostrings::readDNAStringSet(
            paste0("processReads/tzara/", s, "/ASV_", r, ".fasta")
        )
        tzara_regions[[s]][[r]] <- tibble::enframe(
            as.character(tzara_regions[[s]][[r]]),
            name = "seq_id",
            value = r
        )
    }
    tzara_regions[[s]] <-
        purrr::reduce(tzara_regions[[s]], dplyr::full_join, by = "seq_id")
    tzara_regions[[s]] <- dplyr::transmute(
        tzara_regions[[s]],
        seq_id = seq_id,
        `5_8S` = `5_8S`,
        `5_8S_hash` = tzara::seqhash(`5_8S`),
        LSU = stringr::str_c(LSU1, V2, LSU2, V3, LSU3, V4, LSU4),
        LSU_hash = tzara::seqhash(LSU)
    )
    tzara_regions[[s]] <- dplyr::mutate_at(
        tzara_regions[[s]],
        dplyr::vars(dplyr::ends_with("_hash")),
        tidyr::replace_na,
        replace = "missing"
    )
}

#### Combined table ####

# get all the unique 5.8S and LSU sequences
all_5_8S <- dplyr::bind_rows(
    dplyr::select(ampliseq_regions, `5_8S_hash`, `5_8S`),
    purrr::map_dfr(tzara_regions, dplyr::select, `5_8S_hash`, `5_8S`)
)
all_5_8S <- unique(all_5_8S)
all_5_8S <- all_5_8S[complete.cases(all_5_8S),]
all_5_8S <- Biostrings::RNAStringSet(chartr("T", "U", tibble::deframe(all_5_8S)))

all_LSU <- dplyr::bind_rows(
    dplyr::select(ampliseq_regions, LSU_hash, LSU),
    purrr::map_dfr(tzara_regions, dplyr::select, LSU_hash, LSU)
)
all_LSU <- unique(all_LSU)
all_LSU <- all_LSU[complete.cases(all_LSU),]
all_LSU <- Biostrings::RNAStringSet(chartr("T", "U", tibble::deframe(all_LSU)))

#### Alignment ####
# align them independently
# these are pretty quick alignments (~30 min)
# they are a bit faster because we are only aligning the unique sequences
align_5_8S <- DECIPHER::AlignSeqs(all_5_8S, processors = 8)
align_LSU <- DECIPHER::AlignSeqs(all_LSU, processors = 8)

# add an extra sequence to each one for "missing"
# it is just gaps, with the same length as the alignment
align_5_8S["missing"] <-
    Biostrings::RNAStringSet(strrep("-", Biostrings::width(align_5_8S[1])))
align_LSU["missing"] <-
    Biostrings::RNAStringSet(strrep("-", Biostrings::width(align_LSU[1])))

#### Concatenation ####
# now find all unique pairs of 5.8S and LSU
all_both <- dplyr::bind_rows(
    dplyr::select(ampliseq_regions, `5_8S_hash`, LSU_hash),
    purrr::map_dfr(tzara_regions, dplyr::select, `5_8S_hash`, LSU_hash)
)
all_both <- unique(all_both)
# paste together the 5.8S and LSU sequences for each.
align_both <- paste(align_5_8S[all_both$`5_8S_hash`], align_LSU[all_both$LSU_hash], sep = "")
names(align_both) <- paste(all_both$`5_8S_hash`, all_both$LSU_hash, sep = "_")
# make it DNA instead of RNA for fastree
align_both <- Biostrings::DNAStringSet(chartr("U", "T", align_both))

# output the alignment
comparedir <- file.path("processReads", "compare")
if (!dir.exists(comparedir)) dir.create(comparedir)
Biostrings::writeXStringSet(align_both, file.path(comparedir, "comparealn.fasta"))

# remove columns with at least 90% gaps
align_degap <- Biostrings::DNAMultipleAlignment(align_both)
align_degap <- Biostrings::maskGaps(align_degap, min.fraction = 0.9)
align_degap <- as(align_degap, "DNAStringSet")
Biostrings::writeXStringSet(align_degap, file.path(comparedir, "comparealn.degap.fasta"))

#### ML Tree ####
# make a tree with fasttree
# takes about 5 min
system2(
    command = "fasttree",
    args = c("-nt", "-gtr", file.path(comparedir, "comparealn.degap.fasta")),
    stdout = file.path(comparedir, "comparealn.degap.tree")
)

tree <- treeio::read.newick(file.path(comparedir, "comparealn.degap.tree"))
tree <- phangorn::midpoint(tree)

#### Community matrix ####
# make a "community matrix" for the different pipelines
library(magrittr)
ampliseq_reads <- readr::read_tsv(
    "processReads/ampliseq/feature-table.tsv",
    skip = 1,
    col_types = readr::cols(
        .default = readr::col_integer(),
        `#ASV_ID` = readr::col_character()
    )
)

ampliseq_reads <-
    tibble::column_to_rownames(ampliseq_reads, "#ASV_ID") %>%
    rowSums() %>%
    tibble::enframe(name = "seq_id", value = "arrow_ampliseq") %>%
    dplyr::left_join(ampliseq_regions, by = "seq_id") %>%
    dplyr::mutate_at(c("5_8S_hash", "LSU_hash"), tidyr::replace_na, "missing") %>%
    dplyr::transmute(
        seq_id = paste(`5_8S_hash`, LSU_hash, sep = "_"),
        arrow_ampliseq = arrow_ampliseq/sum(arrow_ampliseq)
    ) %>%
    dplyr::group_by(seq_id) %>%
    dplyr::summarize_all(sum)

tzara_reads <- list()
for (s in tzara_sets) {
    s_sym <- rlang::sym(s)
    table <- readRDS(paste0("processReads/tzara/", s, "/asv_table.rds"))
    table <- colSums(table)
    table <- tibble::enframe(table, name = "seq_id", value = s)
    table <- dplyr::left_join(table, tzara_regions[[s]], by = "seq_id")
    table <- dplyr::transmute(
        table,
        seq_id = paste(`5_8S_hash`, LSU_hash, sep = "_"),
        !!s_sym := !!s_sym / sum(!!s_sym)
    )
    table <- dplyr::group_by(table, seq_id)
    table <- dplyr::summarize_all(table, sum)

    tzara_reads[[s]] <- table
}
otu_tab <- purrr::reduce(c(tzara_reads, list(ampliseq_reads)), dplyr::full_join, by = "seq_id")
otu_tab <- dplyr::mutate_if(otu_tab, is.numeric, tidyr::replace_na, 0)
otu_tab <- tibble::column_to_rownames(otu_tab, "seq_id")

otu_tab <- phyloseq::otu_table(otu_tab, taxa_are_rows = TRUE)

#### phyloseq object and UniFrac distances ####
physeq <- phyloseq::phyloseq(tree, otu_tab)
phyloseq::UniFrac(physeq)
phyloseq::UniFrac(physeq, weighted = TRUE)

#### Figure ####
(ggtree::ggtree(tree) +
    ggtree::theme_tree2() +
    ggtree::geom_tiplab(label = "", align = TRUE)) %>%
    ggtree::gheatmap(log(as(otu_tab, "matrix")),
                     colnames_angle = 90,
                     hjust = 1,
                     width = 0.2,
                     legend_title = "log10(abundance)") %>%
    ggplot2::ggsave("processReads/compare/treemap.pdf", plot = ., device = "pdf",
                    width = 8, height = 300, limitsize = FALSE)

physeq_glom <- phyloseq::tip_glom(physeq, h = 0.01)
(ggtree::ggtree(physeq_glom) +
        ggtree::theme_tree2() +
        ggtree::geom_tiplab(label = "", align = TRUE)) %>%
    ggtree::gheatmap(log(as(phyloseq::otu_table(physeq_glom), "matrix")),
                     colnames_angle = 90,
                     hjust = 1,
                     width = 0.2,
                     legend_title = "log10(abundance)") %>%
    ggplot2::ggsave("processReads/compare/treemap_glom.pdf", plot = ., device = "pdf",
                    width = 8, height = 300, limitsize = FALSE)

colSums(phyloseq::otu_table(physeq_glom) > 0)
colSums(phyloseq::otu_table(physeq_glom) > 1/30000) # no global singletons
