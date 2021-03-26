library(targets)
library(tarchetypes)

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
