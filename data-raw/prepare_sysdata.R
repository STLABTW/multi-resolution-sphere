# Generates R/sysdata.rda from the original integral_table2.rds file.
# Run from the package root: Rscript data-raw/prepare_sysdata.R

tab <- readRDS("integral_table2.rds")
stopifnot(is.list(tab), all(c("aa_grid", "tab_vals") %in% names(tab)))

.integral_table <- list(
    aa_grid  = as.numeric(tab$aa_grid),
    tab_vals = as.numeric(tab$tab_vals)
)

save(
    .integral_table,
    file        = "R/sysdata.rda",
    compress    = "xz",
    compression_level = 9
)

message(sprintf(
    "Wrote R/sysdata.rda (length %d, %d KB on disk)",
    length(.integral_table$aa_grid),
    round(file.info("R/sysdata.rda")$size / 1024)
))
