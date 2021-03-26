library(targets)
library(tarchetypes)

outdir <- "output"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
figdir <- file.path(outdir, "figures")
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)
datadir <- file.path(outdir, "data")
if (!dir.exists(datadir)) dir.create(datadir, recursive = TRUE)

comparedir <- file.path("processReads", "compare")
if (!dir.exists(comparedir)) dir.create(comparedir)

# defines lsux_plan
source(here::here("scripts", "01_functions.R"))
# source(here::here("scripts", "02_lsux.R"))
# source(here::here("scripts", "03_dada2.R"))
source(here::here("scripts", "04_otu_tables.R"))
source(here::here("scripts", "05_taxonomy.R"))
source(here::here("scripts", "06_compare_tree.R"))
list(
  otu_table_plan,
  compare_tree_targets,
  taxonomy_targets
)
