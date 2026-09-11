# Genome-wide parallelism in divergence and diversity across three Hydrophis sea snake species

Code, data and figures for the manuscript

> Genome-wide parallelism in divergence and diversity across three Hydrophis sea snake species

Data processing scripts/pipelines are available in the [scripts][scripts] directory.
These are the script responsible for sequence data QC, alignment, variant calling
etc...

Downstream processing scripts are found in the [results][results] directory. Each
of the major analyses have their own directory with a `README.md` in it, which
represents the output of the compiled `.qmd` file. All main text figures/tables
and supplementary figures/tables come from these documents.

- [00-chr2-clean][00-chr2]: Manual curation of chromosome 2 sequence
- [01-alignment-variants][01-aln-var]: Alignment and variant calling summary
- [02-demography][02-dem]: SFS & MSMC2 results
- [03-diversity][03-div]: Within and between species divergence & divergence profiling
- [04-treeness-introgression][04-tree]: Analyses exploring treeness and introgression
- [05-v2r][05-v2r]: V2R content and divergence signals across the chromosome 2 region

The additional directory [00-other][00-other] simply has additional data files
that are used in the scripts above.

Please cite the paper if you find this directory helpful/use any of the code

> TBD

[scripts]: https://github.com/a-lud/hydrophis-popgen/tree/main/scripts
[results]: https://github.com/a-lud/hydrophis-popgen/tree/main/results
[00-chr2]: https://github.com/a-lud/hydrophis-popgen/tree/main/results/00-chr2-clean
[00-other]: https://github.com/a-lud/hydrophis-popgen/tree/main/results/00-other
[01-aln-var]: https://github.com/a-lud/hydrophis-popgen/tree/main/results/01-alignment-variants
[02-dem]: https://github.com/a-lud/hydrophis-popgen/tree/main/results/02-demography
[03-div]: https://github.com/a-lud/hydrophis-popgen/tree/main/results/03-diversity
[04-tree]: https://github.com/a-lud/hydrophis-popgen/tree/main/results/04-treeness-introgression
[05-v2r]: https://github.com/a-lud/hydrophis-popgen/tree/main/results/05-v2r
