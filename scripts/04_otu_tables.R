library(targets)
library(tarchetypes)
library(magrittr)

tar_plan(

    #### Single linkage (Swarm/GeFaST) ####

    tar_file(sl_seq_file, "process/pb_363.swarm.cons.fasta"),
    sl_seqs = load_cons_seqs(sl_seq_file, "swarm_"),

    # load the single-linkage table file and format it
    tar_file(sl_table_file, "process/pb_363_ccs.swarm.table"),
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
        sl_table_xlsx,
        file.path(datadir, "swarm_table.xlsx") %T>%
        openxlsx::write.xlsx(sl_table, .)
    ),

    tar_file(
        sl_table_rds,
        file.path(datadir, "swarm_table.rds") %T>%
        saveRDS(sl_table, .)
    ),

    #### VSEARCH ####

    # load the VSEARCH cluster consensus sequences
    tar_file(vs_seqs_file, "process/pb_363.vclust.cons.fasta"),
    vs_seqs = load_cons_seqs(vs_seqs_file, "OTU_"),

    # load the VSEARCH clustered table file and format it
    tar_file(vs_rawtable_file, "process/pb_363.ccs.vclust.otu_table.txt"),
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
        vs_table_xlsx,
        file.path(datadir, "vsearch_otu_table.xlsx") %T>%
            openxlsx::write.xlsx(vs_table, .)
    ),

    tar_file(
        vs_table_rds,
        file.path(datadir, "vsearch_otu_table.rds") %T>%
            saveRDS(vs_table, .)
    ),

    #### Ampliseq ####
    tar_file(ampliseq_rawtable_file, "processReads/ampliseq/feature-table.tsv"),
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
        ampliseq_table_xlsx,
        file.path(datadir, "ampliseq_table.xlsx") %T>%
        openxlsx::write.xlsx(ampliseq_table, .)
    ),
    tar_file(
        ampliseq_table_rds,
        file.path(datadir, "ampliseq_table.rds") %T>%
            saveRDS(ampliseq_table, .)
    ),

    #### log file
    tar_file(
        otu_log,
        "logs/otu_tables.log" %T>%
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


