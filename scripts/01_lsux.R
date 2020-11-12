# run LSUx on the demultiplexed, primer-free reads, and save the positions.

demux_files <- list.files(
    path = here::here("processReads/demultiplex/demult/soilSamplesFastq"),
    pattern = ".*fastq",
    full.names = TRUE
)

# get the path for the CM which is truncated at the LR5 primer site
# (included in LSUx)
cm_32S_trunc <- system.file(
    file.path("extdata", "fungi_32S_LR5.cm"),
    package = "LSUx"
)

# run on 8 CPUs
start <- Sys.time()
positions <-
    lapply(
        demux_files,
        LSUx::lsux,
        # use the truncated cm which stops at LR5
        cm_32 = cm_32S_trunc,
        # include ITS1
        ITS1 = TRUE,
        cpu = 8,
        # allow 2 Gb ram (per process)
        mxsize = 2048
    )
stop <- Sys.time()
names(positions) <- demux_files

tzara_dir <- here::here("processReads/tzara")
if (!dir.exists(tzara_dir)) dir.create(tzara_dir, recursive = TRUE)
saveRDS(positions, file.path(tzara_dir, "positions.rds"))
