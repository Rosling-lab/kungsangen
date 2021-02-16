# Generate species accumulation curves
library(magrittr)
library(dplyr)
library(ggplot2)
library(phyloseq)

# generates phyloseq objects carbom and carbom_PSH
#source(here::here("R", "phyloseq", "make_phyloseq_obj.R"))

# read the ASV table
community_matrix <- read.table(
    here::here("processReads", "ampliseq", "ampliseq_table.tsv"),
    row.names = 1,
    header = TRUE
)

# read the metadata
samples_df <- read.table(
    here::here("start_files", "meta.txt"),
    header = TRUE,
    row.names = 1,
    stringsAsFactors = FALSE
) %>%
    tidyr::unite("Condition", Sites, Sample_type, sep = "_", remove = FALSE)

# create a phyloseq object
physeq <- phyloseq(
    otu_table(community_matrix, taxa_are_rows = TRUE),
    sample_data(samples_df)
)

outdir <- here::here("output", "figures")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# function to print numbers greater than 1000 with "k"
formatk <- function(x, ...) {
  ifelse(
    x < 1000,
    format(x, ...),
    paste0(format(x / 1000, ...), "k")
  )
}

### Alpha diversity at each sample ####

# maximum number of reads per sample
max(colSums(community_matrix))
min(colSums(community_matrix))

# make an accumulation curve for reads at each sample
seq_accum_sample <-
  # extrapolate the curves up to 2000 reads, with 100 points on each curve
  iNEXT::iNEXT(community_matrix, endpoint = 2000, knots = 100)

# take the mean of all the samples at each
seq_accum_sample_means <-
  # get the plotting data
  fortify(seq_accum_sample) %>%
  # iNEXT doesn't give the same set of x values for each site, so it's hard to
  # take a reasonable average.  So for each site, make a spline function to
  # interpolate the iNEXT results,and get the values at a particular set of x
  # values.
  dplyr::group_by(site) %>%
  dplyr::summarize(
    interp = list(splinefun(x = x, y = y)),
    x = list(c(1, 1:100 * 20)),
    y = list(interp[[1]](x[[1]]))
  ) %>%
  dplyr::select(site, x, y) %>%
  tidyr::unchop(c(x, y)) %>%
  # add wetness annotations
  dplyr::mutate(Type = samples_df[as.character(site), "Sites"]) %>%
  # get the median at each wetness.
  dplyr::group_by(x, Type) %>%
  dplyr::summarize_at("y", mean, .groups = "none")

# make the plot
# get the plotting data
fortify(seq_accum_sample) %>%
  # add forest type annotations
  dplyr::mutate(Type = samples_df[site, "Type"]) %>%
  # plot it
  ggplot(aes(x = x, y = y)) +
  # lines for each sample
  geom_line(aes(
    group = paste(site, method),
    # darker when interpolated, lighter when extrapolated
    alpha = case_when(method == "interpolated" ~ 0.35, TRUE ~0.10)
  )) +
  # use the alpha values I provided as the literal alpha values
  scale_alpha_identity() +
  # small points to represent the actual observed read depth and richness
  geom_point(data = ~filter(., method == "observed"), alpha = 0.4, size = 0.5) +
  # thich dashed line for the mean of all the accumulation curves
  geom_line(data = seq_accum_sample_means, size = 0.8, linetype = "dashed") +
  # seperate facet for each forest type.
  # They are side-by-side so we can compare them visually.
  facet_wrap(facets = ~Type, nrow = 1) +
  xlab("reads") +
  scale_x_continuous() +
  ylab("Species richness")
ggsave(path = outdir, filename = "Accum1.pdf", device = cairo_pdf,
       width = 6.25, height = 3)
ggsave(path = outdir, filename = "Accum1.png", device = "png",
       width = 6.25, height = 3, dpi = 150)

#### Alpha diversity at each condition ####

# Calculate for 4 conditions; i.e. Wet-F, Dry-F, Wet-NF, Dry-NF

#What is the greatest number of reads at a condition?
merge_samples(physeq, "Condition") %>%
  otu_table() %>%
  rowSums() %>%
  max()

# make an accumulation curve for reads at each site
# i.e., what would happen if we had the same samples,
# but sequenced them more (or less)?
seq_accum_site <-
  merge_samples(physeq, "Condition", fun = first) %>%
  otu_table() %>%
  as("matrix") %>%
  t() %>%
  iNEXT::iNEXT(endpoint = 8000, knots = 100)

# make an accumulation curve for samples at each site
# i.e., what would happen if we had the same sequencing depth per sample,
# but got more (or less) samples per site?
samp_accum_site <-
  community_matrix %>%
  vegan::decostand(method = "pa") %>%
  t() %>%
  split.data.frame(samples_df$Condition) %>%
  lapply(t) %>%
  iNEXT::iNEXT(datatype = "incidence_raw", endpoint = 30)

# make plotting data for both types of accumulation curve at each site
accum_site <- list(seq_accum_site, samp_accum_site) %>%
  # get plotting data for each
  lapply(ggplot2::fortify) %>%
  # add a column "xvar" to identify them, and put them together in one tibble
  list(., xvar = c("reads", "samples")) %>%
  purrr::pmap_dfr(tibble::add_column) %>%
  # duplicate the "observed" points as "interpolated" and "extrapolated"
  # so that there isn't a gap in the lines in the plot.
  # this noticable before because they were very closely spaced.
  dplyr::bind_rows(
    dplyr::filter(., method == "observed") %>%
      dplyr::select(-method) %>%
      tidyr::crossing(method = c("interpolated", "extrapolated"))
  ) %>%
    tidyr::separate(
        site,
        into = c("Moisture", "Fritillaria"),
        extra = "merge",
        remove = FALSE
    ) %>%
  # we will label the observed points. but nothing else
  dplyr::mutate(label = ifelse(method == "observed", as.character(Moisture), ""))

# Find the asymptotic species richness estimates
asymp_site <- list(seq_accum_site, samp_accum_site) %>%
  purrr::map(
    ~ cbind(
      dplyr::filter(.$AsyEst, Diversity == "Species richness"),
      .$DataInfo
    ) %>%
      dplyr::select(Estimator, site)
  ) %>%
  # add a column "xvar" to identify them, and put them together in one tibble
  list(., xvar = c("reads", "samples")) %>%
    purrr::pmap_dfr(tibble::add_column) %>%
    tidyr::separate(
        site,
        into = c("Moisture", "Fritillaria"),
        extra = "merge",
        remove = FALSE
    )

## stacked plot

fig_seq_accum_site <-
  dplyr::filter(accum_site, xvar == "reads") %>%
  ggplot(
    aes(x = x, y = y, ymin = y.lwr, ymax = y.upr,
        # within each facet, each site gets its own line.
        group = paste(site, method), color = Moisture, fill = Moisture,
        label = label)
  ) +
  # Dotted horizontal line for the asymptotic estimates
  geom_hline(
    aes(yintercept = Estimator, color = Moisture),
    linetype = "dashed",
    data = dplyr::filter(asymp_site, xvar == "reads"),
    alpha = 0.5
  ) +
  # ribbon for the confidence band
  geom_ribbon(
    alpha = 0.3, # make it semitransparent
    color = NA # don't draw the edges
  ) +
  # line for the estimated curve
  geom_line(aes(alpha = case_when(method == "interpolated" ~ 1, TRUE ~0.5))) +
  scale_alpha_identity() +
  # points for the observed values
  geom_point(data = filter(accum_site, method == "observed", xvar == "reads"),
             alpha = 1) +
    # label the points
    # This is much easier to read than a legend
    ggrepel::geom_text_repel(segment.alpha = 0.2, max.overlaps = 100) +
      facet_wrap( ~ Fritillaria) +
  xlab("reads") +
  scale_x_continuous(breaks = 0:5 * 2000,
                     labels = c(0, paste0(1:5 * 2, "k"))) +
      ylab("Species richness") +
      scale_color_brewer(type = "qual", aesthetics = c("color", "fill"),
                         palette = 2, guide = "none", direction = -1)

fig_samp_accum_site <-
  dplyr::filter(accum_site, xvar == "samples") %>%
  ggplot(
    aes(x = x, y = y, ymin = y.lwr, ymax = y.upr,
        # within each facet, each site gets its own line.
        group = paste(site, method), color = Moisture, fill = Moisture,
        label = label)
  ) +
  # Dotted horizontal line for the asymptotic estimates
  geom_hline(
    aes(yintercept = Estimator, color = Moisture),
    linetype = "dashed",
    data = dplyr::filter(asymp_site, xvar == "samples"),
    alpha = 0.5
  ) +
  # ribbon for the confidence band
  geom_ribbon(
    alpha = 0.3, # make it semitransparent
    color = NA # don't draw the edges
  ) +
  # line for the estimated curve
  geom_line(aes(alpha = case_when(method == "interpolated" ~ 1, TRUE ~0.5))) +
  scale_alpha_identity() +
  # points for the observed values
  geom_point(data = filter(accum_site, method == "observed", xvar == "samples"),
             alpha = 1) +
    # label the points
    # This is much easier to read than a legend
    ggrepel::geom_text_repel(force = 50,
                             segment.alpha = 0.2,
                             seed = 1) +
  # facet by (forest) Type and "xvar"
  # we then move the facet label (aka "strip") for "xvar" to put it at the
  # bottom instead of the top, outside the axis instead of inside it,
  # and remove the shaded box around it.
  # We then allow the different columns to have different ranges on the x axis,
  # and remove the "real" x axis labels.
  # The facet labels have become different x-axis labels for the two columns.
  facet_wrap(~Fritillaria, nrow = 1) +
  xlab("samples") +
  ylab("Species richness") +
  # color brewer has pallettes which are easier to distinguish, although it
  # will always be hard with 8
  scale_color_brewer(type = "qual", aesthetics = c("color", "fill"),
                     palette = 2, guide = "none", direction = -1)

ggpubr::ggarrange(fig_seq_accum_site, fig_samp_accum_site, ncol = 1,
                  labels = "auto")
ggsave(path = outdir, filename = "Accum2.pdf", device = cairo_pdf,
       width = 6.25, height = 6)
ggsave(path = outdir, filename = "Accum2.png", device = "png",
       width = 6.25, height = 6, dpi = 150)
