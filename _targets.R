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

# define which plot types to generate
plot_type_meta <- tibble::tibble(
  ext = c("pdf", "png", "eps"),
  fun = list(rlang::sym("cairo_pdf"), "png", "eps")
)

# defines lsux_plan
source(here::here("scripts", "01_functions.R"))
# source(here::here("scripts", "02_lsux.R"))
# source(here::here("scripts", "03_dada2.R"))
source(here::here("scripts", "04_otu_tables.R"))
source(here::here("scripts", "05_lsux.R"))
source(here::here("scripts", "06_postcluster.R"))
source(here::here("scripts", "07_taxonomy.R"))
source(here::here("scripts", "08_tree.R"))
source(here::here("scripts", "09_accumulation.R"))
source(here::here("scripts", "10_kingdom_tables.R"))
source(here::here("scripts", "11_detection_limits.R"))
source(here::here("scripts", "12_taxonomy_plots.R"))
source(here::here("scripts", "13_venn.R"))
source(here::here("scripts", "14_histograms.R"))
c(
  otu_table_plan,
  lsux_plan,
  its2_cluster_plan,
  taxonomy_plan,
  tree_plan,
  accumulation_plan,
  kingdom_table_plan,
  detection_plan,
  phylotax_plan,
  taxdata_plan,
  taxplot_plan,
  venn_plan,
  hist_plan
)
