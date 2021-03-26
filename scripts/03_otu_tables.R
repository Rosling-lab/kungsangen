library(targets)
library(tarchetypes)
library(magrittr)

tar_plan(

    #### Single linkage (Swarm/GeFaST) ####

    tar_file(sl_seq_file, here::here("process/pb_363.swarm.cons.fasta")),
    sl_seqs = load_sl_seqs(sl_seq_file),

    # load the single-linkage table file and format it
    tar_file(sl_table_file, here::here("process", "pb_363_ccs.swarm.table")),
    sl_rawtable = load_sl_rawtable(sl_table_file, sl_seqs),
    sl_seqtab =
        dplyr::select(sl_rawtable, -OTU) %>%
        tibble::column_to_rownames("seq") %>%
        as.matrix() %>%
        t(),

    sl_nochim = dada2::removeBimeraDenovo(sl_seqtab, verbose = TRUE),

    sl_nosingle = sl_nochim[,colSums(sl_nochim) > 1L],

    sl_table = t(sl_nosingle) %>%
        tibble::as_tibble(rownames = "seq") %>%
        dplyr::left_join(dplyr::select(sl_rawtable, "seq", "OTU"), by = "seq") %>%
        dplyr::select(-seq) %>%
        dplyr::select(OTU, dplyr::everything()) %>%
        dplyr::arrange(readr::parse_number(OTU)),

    tar_file(
        sl_table_tsv,
        here::here("processReads", "swarm_table.tsv") %T>%
        write.table(sl_table, ., sep = "\t")
    ),

    tar_file(
        sl_table_rds,
        here::here("processReads", "swarm_table.rds") %T>%
        saveRDS(sl_table, .)
    ),

    #### VSEARCH ####

    # load the VSEARCH cluster consensus sequences
    tar_file(vs_seqs_file, here::here("processReads", "clusterOTUs", "cluster",
                                      "otus_all20samp.fasta")),
    vs_seqs = Biostrings::readDNAStringSet(vs_seqs_file) %>%
        set_names(., gsub("centroid=(OTU_[0-9]+);.+", "\\1", names(.))) %>%
      as.character(),

    # load the VSEARCH clustered table file and format it
    tar_file(vs_rawtable_file, here::here("processReads", "clusterOTUs",
                                          "cluster", "otu_table.txt")),
    vs_rawtable = readr::read_delim(
        vs_rawtable_file,
        delim = "\t",
        col_types = readr::cols(
            .default = readr::col_integer(),
            `#OTU ID` = readr::col_character()
        )
    ) %>%
        dplyr::left_join(
            tibble::enframe(vs_seqs, name = "#OTU ID", value = "seq")
        ) %>%
        dplyr::group_by(seq) %>%
        dplyr::summarize(
            `#OTU ID` = dplyr::first(`#OTU ID`),
            dplyr::across(where(is.integer), sum)
        ),

    vs_seqtable = dplyr::select(vs_rawtable, -"#OTU ID") %>%
        tibble::column_to_rownames("seq") %>%
        as.matrix(),

    vs_nochim = dada2::removeBimeraDenovo(
        t(vs_seqtable),
        multithread = TRUE,
        verbose = TRUE
    ),

    vs_nosingle = vs_nochim[,colSums(vs_nochim) > 1L],

    vs_table = t(vs_nosingle) %>%
        tibble::as_tibble(rownames = "seq") %>%
        dplyr::left_join(dplyr::select(vs_rawtable, "seq", "#OTU ID"), by = "seq") %>%
        dplyr::select(-seq) %>%
        dplyr::select(OTU = "#OTU ID", dplyr::everything()) %>%
        dplyr::arrange(readr::parse_number(OTU)),

    tar_file(
        vs_table_tsv,
        here::here("processReads/clusterOTUs/vsearch_otu_table.tsv") %T>%
            readr::write_tsv(vs_table, .)
    ),

    tar_file(
        vs_table_rds,
        here::here("processReads/clusterOTUs/vsearch_otu_table.rds") %T>%
            saveRDS(vs_table, .)
    ),

    #### Ampliseq ####
    tar_file(ampliseq_rawtable_file,
             here::here("processReads/ampliseq/feature-table.tsv")),
    ampliseq_table = readr::read_tsv(
        ampliseq_rawtable_file,
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
        dplyr::rename(OTU = "#ASV_ID"),

    tar_file(
        ampliseq_table_tsv,
        here::here("processReads/ampliseq/ampliseq_table.tsv") %T>%
        readr::write_tsv(ampliseq_table, .)
    ),
    tar_file(
        ampliseq_table_rds,
        here::here("processReads/ampliseq/ampliseq_table.rds") %T>%
            saveRDS(ampliseq_table, .)
    ),

    #### log file
    tar_file(
        otu_log,
        here::here("logs/otu_tables.log") %T>%
        writeLines(
            c(
                paste("Single linkage clusters:", sum(sl_seqtab), "reads in",
                      ncol(sl_seqtab), "OTUs"),
                paste("Single linkage clusters after chimera removal:",
                      sum(sl_nochim), "reads in", ncol(sl_nochim), "OTUs"),
                paste("Single linkage clusters after singleton removal:",
                      sum(sl_nosingle), "reads in", ncol(sl_nosingle), "OTUs"),

                paste("VSEARCH 99% clusters:",
                      sum(dplyr::select_if(vs_rawtable, is.integer)),
                      "reads in", nrow(vs_rawtable), "OTUs"),
                paste("VSEARCH 99% clusters after chimera removal:",
                      sum(vs_nochim), "reads in", ncol(vs_nochim), "OTUs"),
                paste("VSEARCH 99% clusters after singleton removal:",
                      sum(vs_nosingle), "reads in", ncol(vs_nosingle), "OTUs")
            ),
            .
        )
    )
) -> otu_table_plan


