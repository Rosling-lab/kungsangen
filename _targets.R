library(targets)
library(tarchetypes)

# defines lsux_plan
# source(here::here("scripts", "01_lsux.R"))
# source(here::here("scripts", "02_dada2.R"))
source(here::here("scripts", "03_otu_tables.R"))
source(here::here("scripts", "04_compare_tree.R"))
source(here::here("scripts", "06_taxonomy.R"))
list(
  otu_table_plan,
  compare_tree_targets,
  taxonomy_targets
)
