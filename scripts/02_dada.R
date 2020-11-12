# use the positions found by LSUx to extract regions from the raw reads and
# run dada2 on each region independently.

if (!exists("ncpu")) ncpu <- 8 # change if using more or less

# load positions from LSUx
# calculated by 01_lsux.R
tzara_dir <- here::here("processReads/tzara")
positions <- readRDS(file.path(tzara_dir, "positions.rds"))

regions <- unique(positions[[1]]$region)
regions

# log to a file
log <- here::here("logs", "tzara.log")
unlink(log)
log <- file(log, "w")
sink(log, type = "output", append = TRUE)
sink(log, type = "message", append = TRUE)

#### Extract regions ####
# based on the positions calculated by LSUx.
# This is not multithreaded, but it only takes a few minutes.
regions_dir <- file.path(tzara_dir, "regions")
for (r in regions) {
    cat("Extracting", r, "\n")
    # files from each region go in their own directory
    the_region_dir <- file.path(regions_dir, r)
    if (!dir.exists(the_region_dir))
        dir.create(the_region_dir, recursive = TRUE)
    # the names of the files are the names of the positions list
    for (f in names(positions)) {
        file_out <- file.path(the_region_dir, basename(f))
        cat("-", f, "->", file_out, "\n")
        # gzip them to take less space
        if (!endsWith(file_out, ".gz")) file_out <- paste0(file_out, ".gz")
        # do the extractions
        tzara::extract_region(
            seq = f,
            positions = positions[[f]],
            region = r,
            outfile = file_out
        )
    }
}

#### Quality filter the regions ####
# Based on mininum and maximum length and maximum expected errors
# Length limits will depend on the region

maxLen <- c(
    ITS1   = 500,
    `5_8S` = 175,
    ITS2   = 500,
    LSU1   = 120,
    V2     = 200,
    LSU2   = 175,
    V3     = 1000,
    LSU3   = 60,
    V4     = 500,
    LSU4   = 175
)

minLen <- c(
    ITS1   = 50,
    `5_8S` = 125,
    ITS2   = 50,
    LSU1   = 90,
    V2     = 125,
    LSU2   = 150,
    V3     = 50,
    LSU3   = 30,
    V4     = 50,
    LSU4   = 150
)

# Use DADA2's filtering functions
# it works on a whole directory of files at a time
filter_dir <- file.path(tzara_dir, "filter")
if (!dir.exists(filter_dir)) dir.create(filter_dir)
for (r in regions) {
    cat("Filtering", r, "\n")
    dada2::filterAndTrim(
        file.path(regions_dir, r),
        file.path(filter_dir, r),
        truncQ = 0,
        minLen = minLen[r],
        maxLen = maxLen[r],
        maxN = 0,
        maxEE = 3,
        multithread = ncpu,
        qualityType = "FastqQuality",
        compress = TRUE,
        verbose = TRUE
    )
}

#### Dereplicate ####
# using dada2
derep <- list()
for (r in regions) {
    cat("Dereplicating", r, "\n")
    derep[[r]] <- dada2::derepFastq(
        file.path(filter_dir, r),
        qualityType = "FastqQuality",
        verbose = TRUE
    )
}

#### Fit error model ####
cat("Fitting error model\n")
err <- dada2::learnErrors(
    fls = derep[["5_8S"]],
    nbases = 1e9,
    errorEstimationFunction = dada2::PacBioErrfun,
    multithread = ncpu,
    verbose = TRUE,
    HOMOPOLYMER_GAP_PENALTY = 0,
    BAND_SIZE = 32,
    pool = TRUE
)

#### run dada on each region ####
dada <- list()
seqtab <- list()
for (r in regions) {
    cat("Denoising", r, "\n")
    dada[[r]] <- dada2::dada(
        derep = derep[[r]],
        err = err,
        pool = TRUE,
        multithread = ncpu,
        verbose = TRUE,
        HOMOPOLYMER_GAP_PENALTY = 0,
        BAND_SIZE = 32
    )
    seqtab[[r]] <- dada2::makeSequenceTable(dada[[r]])
    cat(ncol(seqtab[[r]]), "ASVs in", sum(seqtab[[r]]), "sequences\n")
}

#### generate dadamaps ####
dadamaps <- list()
for (r in regions) {
    dadamaps[[r]] <- tzara::dadamap(
        derep = derep[[r]],
        dada = dada[[r]],
        filename = file.path(filter_dir, r, names(derep[[r]]))
    )
}

#### get denoised regions for each sequence ####
allseqs <- tzara::reconstruct(dadamaps)

#### dereplicate the denoised individual reads ####
allseqs <- dplyr::group_by_at(allseqs, dplyr::vars(ITS1:LSU4))
allseqs <- dplyr::summarize(allseqs, nread = dplyr::n())
#### get the consensus for non-ITS2 regions ####
library(magrittr)
keycol <- "ITS2"
othercols <- setdiff(regions, keycol)

preconseqs <-
    tidyr::pivot_longer(
        allseqs,
        cols = all_of(othercols),
        names_to = "region",
        values_to = "seq"
    ) %>%
    dplyr::filter(!is.na(.data[[keycol]])) %>%
    dplyr::group_by_at(c(keycol, "region", "seq")) %>%
    dplyr::summarize_at("nread", sum) %>%
    dplyr::group_by_at(c(keycol, "region"))

conseqs <-
    preconseqs %>%
    dplyr::summarize(
        seq = tzara::cluster_consensus(
            seq,
            nread = nread,
            ncpus = ncpu,
            simplify = TRUE
        ),
        nread = sum(nread)
    ) %>%
    tidyr::pivot_wider(names_from = "region", values_from = "seq") %>%
    dplyr::ungroup() %>%
    dplyr::filter(complete.cases(.)) %>%
    dplyr::arrange(dplyr::desc(nread)) %>%
    dplyr::mutate(
        full = do.call(paste0, .[regions]),
        name = paste0("ASV_", formatC(seq_along(full), width = 4, flag = "0"))
    )

consensus <- list()
for (r in c(regions, "full")) {
    consensus[[r]] <- tibble::deframe(conseqs[c("name", r)]) %>%
        chartr("U", "T", .) %>%
        Biostrings::DNAStringSet()

    Biostrings::writeXStringSet(
        consensus[[r]],
        file.path(tzara_dir, paste0("ASV_", r, ".fasta"))
    )
}

asvtab <- seqtab$ITS2[,as.character(consensus$ITS2)]
cat("Calculated", ncol(asvtab), "consensus ASVs covering", sum(asvtab),
    "reads.\n")
colnames(asvtab) <- names(consensus$ITS2)
rownames(asvtab) <- sub(".fastq.gz", "", rownames(asvtab), fixed = TRUE)
saveRDS(asvtab, file.path(tzara_dir, "asv_table.rds"))
readr::write_csv(
    tibble::as_tibble(asvtab, rownames = "sample"),
    file.path(tzara_dir, "asv_table.csv")
)

sink(type = "output")
sink(type = "message")
