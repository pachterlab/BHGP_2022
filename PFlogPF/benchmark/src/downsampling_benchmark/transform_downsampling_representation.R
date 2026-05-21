pa <- argparser::arg_parser("Take a matrix and generate a low-dimensional representation")
pa <- argparser::add_argument(pa, "--transformation", type = "character", help = "The name of the transformation")
pa <- argparser::add_argument(pa, "--input_id", type = "character", help = "The id of a file in output/results")
pa <- argparser::add_argument(pa, "--data_mode", type = "character", help = "Either 'full' or 'reduced'")
pa <- argparser::add_argument(pa, "--pca_dim", type = "numeric", help = "The dimension of the returned representation")
pa <- argparser::add_argument(pa, "--alpha", type = "character", default = "FALSE", help = "The alpha parameter. Ignored by some transformations.")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)
stopifnot(pa$data_mode %in% c("full", "reduced"))

source("src/transformations/transformation_helper.R")

dataset <- readRDS(file.path(pa$working_dir, "results", pa$input_id))
UMI <- dataset[[pa$data_mode]]
expressed_cells <- matrixStats::colSums2(UMI) > 0
expressed_genes <- matrixStats::rowSums2(UMI) > 0
UMI <- UMI[expressed_genes, expressed_cells, drop = FALSE]

alpha <- pa$alpha
if(pa$alpha == "global"){
  alpha <- "global"
}else if(! is.na(suppressWarnings(readr::parse_double(pa$alpha, na = character(0L))))){
  alpha <- readr::parse_double(pa$alpha)
}else if(! is.na(suppressWarnings(readr::parse_logical(pa$alpha, na = character(0L))))){
  alpha <- readr::parse_logical(pa$alpha)
}else{
  stop("Cannot parse alpha=", alpha)
}

sf <- MatrixGenerics::colSums2(UMI)
sf <- sf / mean(sf)

duration <- system.time({
  trans_dat <- all_transformations[[pa$transformation]](UMI, sf, alpha)
  rep <- make_lowdim_representation(pa$transformation, trans_dat, pa$pca_dim)
})

write.table(
  data.frame(name = names(duration), seconds = as.vector(duration)),
  file.path(pa$working_dir, "duration", pa$result_id),
  sep = "\t", row.names = FALSE, quote = FALSE
)
saveRDS(rep, file.path(pa$working_dir, "results", pa$result_id))
