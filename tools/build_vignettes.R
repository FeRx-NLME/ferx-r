root <- normalizePath(".")
doc_dir <- file.path(root, "inst", "doc")
dir.create(doc_dir, recursive = TRUE, showWarnings = FALSE)

vigs <- list.files(file.path(root, "vignettes"), pattern = "\\.Rmd$", full.names = TRUE)
for (f in vigs) {
  out <- file.path(doc_dir, sub("\\.Rmd$", ".html", basename(f)))
  idx <- file.path(doc_dir, sub("\\.Rmd$", ".R",    basename(f)))
  message("Rendering: ", basename(f))
  rmarkdown::render(f, output_file = out, quiet = TRUE)
  knitr::purl(f, output = idx, quiet = TRUE)
}
message("Done")
