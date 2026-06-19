library(tidyverse)

pa <- argparser::arg_parser("Make a table of the siRNA quantitative recovery benchmark results")
pa <- argparser::add_argument(pa, "--sirna_quantitative_recovery_results", type = "character", nargs = Inf, help = "A list of siRNA quantitative recovery result ids")
pa <- argparser::add_argument(pa, "--working_dir", type = "character", help = "The directory that contains the params, results, scripts etc.")
pa <- argparser::add_argument(pa, "--result_id", type = "character", help = "The result_id")
pa <- argparser::parse_args(pa)

print(pa)

res <- lapply(pa$sirna_quantitative_recovery_results, function(res_id){
  res <- read_tsv(file.path(pa$working_dir, "results", res_id), show_col_types = FALSE)
  duration <- map_df(res$representation_result_id, function(id){
    dur <- read_tsv(file.path(pa$working_dir, "duration", id), show_col_types = FALSE)
    tibble(
      cputime_sec = sum(deframe(dur)[c("user.self", "sys.self", "user.child", "sys.child")]),
      elapsed_sec = deframe(dur)["elapsed"]
    )
  })
  bind_cols(res, duration)
}) %>%
  bind_rows()

write_tsv(res, file.path(pa$working_dir, "results", pa$result_id))
