lvls <- paste0("chr", 1:15)

dir_ls("results/03-diversity/within", glob = "*tajima*") |>
    read_tsv(col_types = cols(), id = "window_size") |>
    mutate(
        chromosome = factor(chromosome, levels = lvls),
        window_size = str_extract(window_size, "100kb|50kb")
    ) |>
    arrange(chromosome, window_pos_1) |>
    group_by(pop, window_size) |>
    group_split() |>
    walk(\(x) {
        id <- unique(x$pop)
        w <- unique(x$window_size)

        x |>
            write_tsv(file = glue::glue("results/03-diversity/{id}-{w}_tajima_d.txt"))
    })

dir_ls("results/03-diversity/within", glob = "*pi*") |>
    read_tsv(col_types = cols(), id = "window_size") |>
    mutate(
        chromosome = factor(chromosome, levels = lvls),
        window_size = str_extract(window_size, "100kb|50kb")
    ) |>
    arrange(chromosome, window_pos_1) |>
    group_by(pop, window_size) |>
    group_split() |>
    walk(\(x) {
        id <- unique(x$pop)
        w <- unique(x$window_size)

        x |>
            write_tsv(file = glue::glue("results/03-diversity/{id}-{w}_pi.txt"))
    })

dir_ls("results/03-diversity/within", glob = "*theta*") |>
    read_tsv(col_types = cols(), id = "window_size") |>
    mutate(
        chromosome = factor(chromosome, levels = lvls),
        window_size = str_extract(window_size, "100kb|50kb")
    ) |>
    arrange(chromosome, window_pos_1) |>
    group_by(pop, window_size) |>
    group_split() |>
    walk(\(x) {
        id <- unique(x$pop)
        w <- unique(x$window_size)

        x |>
            write_tsv(file = glue::glue("results/03-diversity/{id}-{w}_watterson_theta.txt"))
    })
