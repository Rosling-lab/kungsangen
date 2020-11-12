# use the positions found by LSUx to extract regions from the raw reads and
# run dada2 on each region independently.

tzara_dir <- here::here("processReads/tzara")
positions <- readRDS(file.path(tzara_dir, "positions.rds"))

regions <- unique(positions[[1]]$region)
regions
