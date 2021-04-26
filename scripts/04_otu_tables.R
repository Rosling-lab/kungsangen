library(targets)
library(tarchetypes)
library(magrittr)

tar_plan(
  tar_map(
    values = list(
      cluster_type = c("swarm", "vclust"),
      cluster_prefix = c("swarm_", "vclust_"),
      id = c("sl", "vs")
    ),
    # load consensus sequences
    tar_file(seq_file, paste0("process/pb_363.", cluster_type, ".cons.fasta")),
    tar_target(seqs, load_cons_seqs(seq_file, cluster_prefix)),
    # load OTU table
    tar_file(table_file, paste0("process/pb_363_ccs.", cluster_type, ".table")),
    tar_target(rawtable, load_rawtable(table_file, seqs)),
    # format OTU table for dada2
    tar_target(
      seqtab,
      dplyr::select(rawtable, -OTU) %>%
        dplyr::group_by(seq) %>%
        dplyr::summarize_all(sum) %>%
        tibble::column_to_rownames("seq") %>%
        as.matrix() %>%
        t()
    ),
    # remove chimeras
    tar_target(
      nochim,
      dada2::removeBimeraDenovo(seqtab, verbose = TRUE, multithread = TRUE)
    ),
    # remove singletons
    tar_target(nosingle, nochim[,colSums(nochim) > 1L]),
    # output OTU table
    tar_target(
      table,
      t(nosingle) %>%
        tibble::as_tibble(rownames = "seq") %>%
        dplyr::left_join(dplyr::select(rawtable, "seq", "OTU"), by = "seq") %>%
        dplyr::select(-seq) %>%
        dplyr::select(OTU, dplyr::everything()) %>%
        dplyr::arrange(readr::parse_number(OTU))
    ),
    tar_map(
      values = list(type = c("rds", "xlsx")),
      tar_file(
        tableout,
        write_and_return_file(
          table,
          file.path(datadir, paste0(cluster_type, "_table.", type)),
          type = type
        )
      )
    ),
    names = id
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
  tar_map(
    values = list(type = c("rds", "xlsx")),
    tar_file(
      ampliseq_tableout,
      write_and_return_file(
        ampliseq_table,
        file.path(datadir, paste0("ampliseq_table", type)),
        type = type
      )
    )
  ),
   #### log file
  tar_file(
      otu_log,
      write_and_return_file(
          c(
              paste("Single linkage clusters:", sum(seqtab_sl), "reads in",
                    ncol(seqtab_sl), "OTUs"),
              paste("Single linkage clusters after chimera removal:",
                    sum(nochim_sl), "reads in", ncol(nochim_sl), "OTUs"),
              paste("Single linkage clusters after singleton removal:",
                    sum(nosingle_sl), "reads in", ncol(nosingle_sl), "OTUs"),
              paste("VSEARCH 99% clusters:",
                    sum(dplyr::select_if(rawtable_vs, is.integer)),
                    "reads in", nrow(rawtable_vs), "OTUs"),
              paste("VSEARCH 99% clusters after chimera removal:",
                    sum(nochim_vs), "reads in", ncol(nochim_vs), "OTUs"),
              paste("VSEARCH 99% clusters after singleton removal:",
                    sum(nosingle_vs), "reads in", ncol(nosingle_vs), "OTUs")
          ),
          "logs/otu_tables.log",
          "txt"
      )
  )
) -> otu_table_plan
