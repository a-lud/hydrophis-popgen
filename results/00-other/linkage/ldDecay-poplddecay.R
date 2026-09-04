# ---------------------------------------------------------------------------- #
# Plot LD-decay from PopLDdecay
#
# This code is an adapted version of the code produced by PopLDdecay. I've
# simply used ggplot2 to generate the figures rather than base-R.

suppressPackageStartupMessages({
    library(tidyverse)
    library(ggtext)
    library(here)
    library(fs)
    library(ragg)
})

# Colour profiles per species
pal_hma <- "#fca311"
pal_hcu <- "#59A14FFF"
pal_hel <- "#4E79A7FF"
lvls <- c("HMA", "HCU", "HEL")

# ---------------------------------------------------------------------------- #
# Data
df.ld <- dir_ls(here("results", "00-other", "linkage"), glob = "*.bin.gz") |>
    map(
        read_tsv,
        col_names = c(
            "dist", "mean_r2", "mean_D",
            "sum_r2", 'sum_D', "nPairs"
        ),
        col_types = cols(),
        skip = 1
    ) |>
    list_rbind(names_to = 'Species') |>
    mutate(
        Species = toupper(str_remove(basename(Species), ".bin.gz")),
        Species = factor(Species, levels = lvls)
    )

# ---------------------------------------------------------------------------- #
# Decay plots
p <- df.ld |>
    ggplot(aes(x = dist, y = mean_r2, colour = Species)) +
    geom_line(linewidth = 1.5) +
    geom_vline(
        xintercept = c(25e3, 100e3),
        linetype = c("dashed", "solid")
    ) +
    scale_colour_manual(values = c(pal_hma, pal_hcu, pal_hel)) +
    scale_x_continuous(
        breaks = scales::breaks_pretty(n = 10),
        labels = scales::label_number(scale = 1e-3)
    ) +
    scale_y_continuous(
        breaks = scales::breaks_pretty(n = 6),
        limits = c(0, 0.6)
    ) +
    labs(
        x = "\nDistance (kb)",
        y = "Mean r<sup>2</sup><br>",
    ) +
    theme_classic() +
    theme(
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 18),
        axis.title.y = element_markdown(size = 18),
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 16),
        legend.position = c(0.9, 0.9)
    )

p

agg_png(
    filename = here(
        "manuscript",
        "figures",
        "supplementary",
        "supp-figure-X-ld-curves.png"
    ),
    width = 1000,
    height = 1000,
    units = "px",
    scaling = 1.2
)
p
invisible(dev.off())
