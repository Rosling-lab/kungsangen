# run LSUx on the demultiplexed, primer-free reads, and save the positions.
# this is only necessary for the tzara ASVs (which we aren't using)
library(targets)
library(tarchetypes)

lsux_plan <-
    tar_plan(
        # find all the demultiplexed fastq files
        tar_files(
            demux_files,
            list.files(
                path = here::here("processReads/demultiplex/demult/soilSamplesFastq"),
                pattern = ".*fastq",
                full.names = TRUE
            )
        ),

        # get the path for the CM which is truncated at the LR5 primer site
        # (included in LSUx)
        tar_file(
            cm_32S_trunc,
            system.file(
                file.path("extdata", "fungi_32S_LR5.cm"),
                package = "LSUx"
            )
        ),

        # run on 1 CPU (but these can be done in parallel)
        tar_fst_tbl(
            positions,
            LSUx::lsux(
                demux_files,
                # use the truncated cm which stops at LR5
                cm_32 = cm_32S_trunc,
                # include ITS1
                ITS1 = TRUE,
                cpu = 8,
                # allow 2 Gb ram (per process)
                mxsize = 2048
            ) %>%
                tibble::add_column(file = demux_files),
            pattern = map(demux_files)
        )
    )
