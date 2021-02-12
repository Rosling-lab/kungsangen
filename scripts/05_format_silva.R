# make dada2 and SINTAX-compatible versions of the SILVA LSU NR99 database.
# the base file is generated from
# SILVA_138.1_LSURef_NR99_tax_silva_trunc.fasta.gz
# using sed '/^>/!y/uU/tT/; /^>/s/^>[^ ] //'

# load just the fasta headers
ref <- Biostrings::readDNAStringSet(
    here::here("reference", "SILVA_138.1_LSURef_NR99_tax_silva_dna_euk_trim.fasta.gz")
)

# parse the headers into accession number and taxonomy
headers <- stringr::str_split_fixed(names(ref), " ", 2)
colnames(headers) <- c("accno", "taxonomy")
headers <- tibble::as_tibble(headers)
headers <- tidyr::extract(
    headers,
    "taxonomy",
    into = c("taxonomy", "species"),
    regex = "(.+);([^;]+)$"
)

# This function taxes a character vector, and returns a vector of the same
# length, cumulatively joining the elements of the vector, with ";" as a separator.
# e.g., c("Fungi", "Basidiomycota", "Agaricomycetes", "Sebacinales") ->
# c("Fungi", "Fungi;Basidiomycetes", "Fungi;Basidiomycetes;Agaricomycetes",
#   "Fungi;Basidiomycetes;Agaricomycetes;Sebacinales")
cumcollapse <- function(x) {
    purrr::map_chr(
        seq_along(x),
        ~paste(x[1:.], collapse = ";")
    )
}

# get all the distinct (terminal) taxa
all_taxa <- dplyr::distinct(headers, taxonomy)

# for each terminal taxon, make a list of all its supertaxa,
# and fully qualify each of them
# don't include "Incertae Sedis" taxa; those will be recreated later.
all_taxa$taxon <- all_taxa$taxonomy
all_taxa <- tidyr::separate_rows(all_taxa, taxon, sep = ";")
all_taxa <- dplyr::filter(all_taxa, !endsWith(taxon, "Incertae Sedis"))
all_taxa <- dplyr::group_by(all_taxa, taxonomy)
all_taxa <- dplyr::mutate(all_taxa, full_taxon = cumcollapse(taxon))
all_taxa <- dplyr::mutate(all_taxa, full_taxon = paste0(full_taxon, ";"))

# load the taxonomy
taxa <- readr::read_delim(
    here::here("reference", "tax_slv_lsu_138.1.txt.gz"),
    delim = "\t",
    col_names = c("taxonomy", "taxid", "rank", "remark", "version"),
    col_types = "ccccc"
)

# all the ranks included in the taxonomy, in order from highest to lowest
silva_ranks <- c("domain", "major_clade", "superkingdom", "kingdom",
                 "subkingdom", "superphylum", "phylum", "subphylum",
                 "infraphylum", "superclass", "class", "subclass", "infraclass",
                 "superorder", "order", "suborder", "superfamily", "family",
                 "subfamily", "genus")

# add rank information to all the taxa
all_taxa <- dplyr::left_join(all_taxa, taxa, by = c("full_taxon" = "taxonomy"))
all_taxa <- dplyr::mutate(
    all_taxa,
    rank = factor(rank, levels = silva_ranks, ordered = TRUE)
)

# make a table giving the ranked classification for each terminal taxon
all_taxa <- dplyr::select(all_taxa, taxonomy, taxon, rank)
all_taxa <- dplyr::group_by(all_taxa, taxonomy, rank)
all_taxa <- dplyr::summarise(all_taxa, taxon = dplyr::first(taxon))
all_taxa <- dplyr::ungroup(all_taxa)
all_taxa <- tidyr::pivot_wider(all_taxa, names_from = "rank", values_from = "taxon")

# make values for all the primary ranks using given values
# when there is no given value, take the name from the next
# higher rank where there is a value.
# If no lower rank is given, use "unidentified_"
# If there is a lower rank given, use "_Incertae_sedis"
unknown = "unidentified_"
incertae = "_Incertae_sedis"
interpolate_rank <- function(x, x2, genus2) {
    dplyr::coalesce(
        x,
        ifelse(
            x2 == genus2,
            paste0(unknown, x2),
            paste0(x2, incertae)
        )
    )
}
all_taxa <- dplyr::mutate(
    all_taxa,
    kingdom2 = dplyr::coalesce(kingdom, superkingdom, major_clade, domain),
    phylum2 = dplyr::coalesce(phylum, superphylum, subkingdom, kingdom2),
    class2 = dplyr::coalesce(class, superclass, infraphylum, subphylum, phylum2),
    order2 = dplyr::coalesce(order, superorder, infraclass, subclass, class2),
    family2 = dplyr::coalesce(family, superfamily, suborder, order2),
    genus2 = dplyr::coalesce(genus, subfamily, family2),
    genus = dplyr::coalesce(genus, paste0(unknown, genus2)),
    family = interpolate_rank(family, family2, genus2),
    order = interpolate_rank(order, order2, genus2),
    class = interpolate_rank(class, class2, genus2),
    phylum = interpolate_rank(phylum, phylum2, genus2),
    kingdom = interpolate_rank(kingdom, kingdom2, genus2),
    dada2 = paste(domain, kingdom, phylum, class, order, family, genus, sep = ";"),
    sintax = glue::glue("tax=d:{domain},k:{kingdom},p:{phylum},c:{class},o:{order},f:{family},g:{genus}")
)

# output in dada2 and sintax formats

headers <- dplyr::left_join(
    headers,
    dplyr::select(all_taxa, taxonomy, dada2, sintax),
    by = "taxonomy"
)

ref_dada <- magrittr::set_names(ref, headers$dada2)
Biostrings::writeXStringSet(
    ref_dada,
    here::here("reference", "silva.LSU.dada2.fasta.gz"),
    compress = TRUE
)
ref_sintax <- magrittr::set_names(ref, headers$sintax)
Biostrings::writeXStringSet(
    ref_sintax,
    here::here("reference", "silva.LSU.sintax.fasta.gz"),
    compress = TRUE
)

# clean up
remove(ref, ref_dada, ref_sintax, headers, taxa, all_taxa)
gc()
