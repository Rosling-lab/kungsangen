# get all unique sequences
# much of code duplicated from 03_compare_tree.R

ampli_seqs <- readr::read_tsv(
    "processReads/ampliseq/qiime2_ASV_table.tsv",
    col_types = readr::cols(
        .default = readr::col_double(),
        Feature.ID = readr::col_character(),
        Taxon = readr::col_character(),
        sequence = readr::col_character()
    )
)

ampli_seqs <- tibble::deframe(ampli_seqs[,c("Feature.ID", "sequence")])

sl_seqs <- Biostrings::readDNAStringSet(
    here::here("process/pb_363.swarm.cons.fasta")
)

sl_seqs <- as.character(sl_seqs)

vs_seqs <- Biostrings::readDNAStringSet(
    here::here("processReads", "clusterOTUs", "cluster", "otus_all20samp.fasta")
)

vs_seqs <- as.character(vs_seqs)

# combine all representative sequences from all three clustering algorithms
all_seqs <- c(ampli_seqs, sl_seqs, vs_seqs)
all_seqs <- unique(all_seqs)
names(all_seqs) <- tzara::seqhash(all_seqs)

# get the path for the CM which is truncated at the LR5 primer site
# (included in LSUx)
cm_32S_trunc <- system.file(
    file.path("extdata", "fungi_32S_LR5.cm"),
    package = "LSUx"
)

# divide up the sequences
positions <- LSUx::lsux(
    all_seqs,
    cm_32S = cm_32S_trunc,
    ITS1 = TRUE,
    cpu = 8,
    # allow 2 Gb ram (per process)
    mxsize = 2048
)

all_regions <- purrr::map2(
    .x = c("ITS1", "LSU1"),
    .y = c("ITS2", "LSU4"),
    tzara::extract_region,
    seq = all_seqs,
    positions = positions
)
names(all_regions) <- c("ITS", "LSU")

tax_silva_full <- phylotax::taxonomy_sintax(
    all_seqs,
    here::here("reference", "silva.LSU.sintax.fasta.gz"),
    min_confidence = 80,
    multithread = 8
)

tax_rdp_full <- phylotax::taxonomy_sintax(
    all_seqs,
    here::here("reference", "rdp_train.LSU.sintax.fasta.gz"),
    min_confidence = 80,
    multithread = 8
)

tax_unite_full <- phylotax::taxonomy_sintax(
    all_seqs,
    here::here("reference", "unite.ITS.sintax.fasta.gz"),
    min_confidence = 80,
    multithread = 8
)

tax_silva_trim <- phylotax::taxonomy_sintax(
    all_regions$LSU,
    here::here("reference", "silva.LSU.sintax.fasta.gz"),
    min_confidence = 80,
    multithread = 8
)

tax_rdp_trim <- phylotax::taxonomy_sintax(
    all_regions$LSU,
    here::here("reference", "rdp_train.LSU.sintax.fasta.gz"),
    min_confidence = 80,
    multithread = 8
)

tax_unite_trim <- phylotax::taxonomy_sintax(
    all_regions$ITS,
    here::here("reference", "unite.ITS.sintax.fasta.gz"),
    min_confidence = 80,
    multithread = 8
)

tax_trim <- dplyr::bind_rows(
    phylotax::taxtable_sintax(tax_unite_trim, min.confidence = 80) %>%
        tibble::add_column(region = "ITS", method = "unite"),
    phylotax::taxtable_sintax(tax_silva_trim, min.confidence = 80) %>%
        tibble::add_column(region = "LSU", method = "silva")
) %>%
    dplyr::filter(confidence >= 0.8)

tax_both <-
    dplyr::filter(
        tax_trim,
        !startsWith(taxon, "unidentified"),
        !endsWith(taxon, "ncertae_sedis")
    ) %>%
    dplyr::group_by(label, rank) %>%
    dplyr::filter(dplyr::n() > 1) %>%
    dplyr::summarize(match = dplyr::n_distinct(taxon) == 1) %>%
    dplyr::ungroup() %>%
    dplyr::add_count(rank) %>%
    dplyr::filter(match) %>%
    dplyr::group_by(rank, n) %>%
    dplyr::summarize(frac = dplyr::n() / unique(n))

cumcollapse <- function(x) {
    purrr::map_chr(
        seq_along(x),
        ~paste(x[1:.], collapse = ";")
    )
}

tax_compare <-
    tax_trim %>%
    dplyr::group_by(label, region) %>%
    dplyr::arrange(rank) %>%
    dplyr::mutate(c12n = cumcollapse(taxon)) %>%
    dplyr::group_by(label, rank) %>%
    dplyr::filter(
        !startsWith(taxon, "unidentified"),
        !endsWith(taxon, "ncertae_sedis")
    ) %>%
    dplyr::filter(dplyr::n_distinct(taxon) > 1) %>%
    dplyr::select(label, taxon, rank, method, c12n) %>%
    tidyr::pivot_wider(names_from = "method", values_from = c("taxon", "c12n")) %>%
    dplyr::ungroup() %>%
    dplyr::count(rank, taxon_unite, taxon_silva, c12n_unite, c12n_silva)
