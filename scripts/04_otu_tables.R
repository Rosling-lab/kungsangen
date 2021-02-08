library(magrittr)

sl_seqs <- Biostrings::readDNAStringSet(
    here::here("process/pb_363.swarm.cons.fasta")
)

names(sl_seqs) <- gsub("consensus", "swarm_", names(sl_seqs))

# load the single-linkage table file and format it
sl_table <- readr::read_delim(
    here::here("process", "pb_363_ccs.swarm.table"),
    delim = " ",
    col_names = c("Sample", "OTU", "reads"),
    col_types = "cci"
) %>%
    dplyr::arrange(OTU) %>%
    tidyr::pivot_wider(names_from = "Sample", values_from = "reads", values_fill = 0L) %>%
    dplyr::left_join(tibble::enframe(as.character(sl_seqs), name = "OTU", value = "seq"), by = "OTU")

sl_seqtab <- sl_table %>%
    dplyr::select(-OTU) %>%
    tibble::column_to_rownames("seq") %>%
    as.matrix() %>%
    t()

cat("Single linkage clusters:", sum(sl_seqtab), "reads in", ncol(sl_seqtab), "OTUs")

sl_nochim <- dada2::removeBimeraDenovo(sl_seqtab, verbose = TRUE)

cat("Single linkage clusters after chimera removal:", sum(sl_nochim), "reads in", ncol(sl_nochim), "OTUs")

sl_nosingle <- sl_nochim[,colSums(sl_nochim) > 1L]

cat("Single linkage clusters after singleton removal:", sum(sl_nosingle), "reads in", ncol(sl_nosingle), "OTUs")

sl_nosingle <- t(sl_nosingle) %>%
    tibble::as_tibble(rownames = "seq") %>%
    dplyr::left_join(dplyr::select(sl_table, "seq", "OTU"), by = "seq") %>%
    dplyr::select(-seq) %>%
    dplyr::select(OTU, dplyr::everything()) %>%
    dplyr::arrange(readr::parse_number(OTU))

write.table(sl_nosingle, here::here("processReads", "swarm_table.tsv"),
            sep = "\t")
saveRDS(sl_nosingle, here::here("processReads", "swarm_table.rds"))


# load the VSEARCH cluster consensus sequences
vs_seqs <- Biostrings::readDNAStringSet(
    here::here("processReads", "clusterOTUs", "cluster", "otus_all20samp.fasta")
)

# load the VSEARCH clustered table file and format it
vs_table <- readr::read_delim(
    here::here("processReads/clusterOTUs/cluster/otu_table.txt"),
    delim = "\t",
    col_types = readr::cols(
        .default = readr::col_integer(),
        `#OTU ID` = readr::col_character()
    )
) %>%
    dplyr::left_join(tibble::enframe(as.character(vs_seqs), name = "#OTU ID", value = "seq")) %>%
    dplyr::group_by(seq) %>%
    dplyr::summarize(
        `#OTU ID` = dplyr::first(`#OTU ID`),
        dplyr::across(where(is.integer), sum)
    )

cat("VSEARCH 99% clusters:", sum(dplyr::select_if(vs_table, is.integer)), "reads in", nrow(vs_table), "OTUs")

vs_seqtable <- dplyr::select(vs_table, -"#OTU ID") %>%
    tibble::column_to_rownames("seq") %>%
    as.matrix()

vs_nochim <- dada2::removeBimeraDenovo(
    t(vs_seqtable),
    multithread = TRUE,
    verbose = TRUE
)

cat("VSEARCH 99% clusters after chimera removal:",
    sum(vs_nochim), "reads in", ncol(vs_nochim), "OTUs")

vs_nosingle <- vs_nochim[,colSums(vs_nochim) > 1L]
cat("VSEARCH 99% clusters after singleton removal:", sum(vs_nosingle), "reads in", ncol(vs_nosingle), "OTUs")

vs_nosingle <- t(vs_nosingle) %>%
    tibble::as_tibble(rownames = "seq") %>%
    dplyr::left_join(dplyr::select(vs_table, "seq", "#OTU ID"), by = "seq") %>%
    dplyr::select(-seq) %>%
    dplyr::select(OTU = "#OTU ID", dplyr::everything()) %>%
    dplyr::arrange(readr::parse_number(OTU))

readr::write_tsv(vs_nosingle,
                 here::here("processReads/clusterOTUs/vsearch_otu_table.tsv"))
saveRDS(vs_nosingle,
        here::here("processReads/clusterOTUs/vsearch_otu_table.rds"))

## Ampliseq

ampliseq_table <- readr::read_tsv(
    here::here("processReads/ampliseq/feature-table.tsv"),
    skip = 1,
    col_types = readr::cols(
        .default = readr::col_integer(),
        `#ASV_ID` = readr::col_character()
    )) %>%
    dplyr::rename_if(
        is.integer,
        ~stringr::str_match(., "([1-5g])([EW])([FN])S.*") %>%
            as.data.frame() %>%
            dplyr::mutate(
                V2 = ifelse(V2 == "g", "5", V2),
                V3 = ifelse(V3 == "E", "Dry", "Wet")
            ) %$%
            paste(V3, V4, V2, sep = "_")
    ) %>%
    dplyr::rename(OTU = "#ASV_ID")

readr::write_tsv(ampliseq_table,
                 here::here("processReads/ampliseq/ampliseq_table.tsv"))
saveRDS(ampliseq_table,
        here::here("processReads/ampliseq/ampliseq_table.rds"))
