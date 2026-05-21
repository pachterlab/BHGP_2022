library(tidyverse)

.OUTPUT_FOLDER <- normalizePath("output", mustWork = FALSE)
source("src/job_management_utils.R")

instructions <- yaml::read_yaml("job_overview.yaml")
di4 <- instructions[["sirna_quantitative_recovery"]]

get_data_for_downsampling_benchmark <- function(dataset_name, seed){
  wrap_script(
    "src/downsampling_benchmark/download_deeply_sequenced_datasets.R",
    params = list(dataset = dataset_name, seed = seed),
    duration = "0-01:00:00", memory = "40GB"
  )
}

transform_deeply_sequenced_representation <- function(transformation, data, data_mode = c("full", "reduced"), pca_dim, alpha){
  data_mode <- match.arg(data_mode)
  params <- list(transformation = transformation, input_id = data$result_id, data_mode = data_mode, pca_dim = pca_dim, alpha = alpha)
  if(transformation %in% c("dino", "glmpca", "newwave")){
    duration <- "2-00:00:00"
    n_cpus <- 1
  }else{
    duration <- "0-05:00:00"
    n_cpus <- 1
  }
  wrap_script(
    "src/downsampling_benchmark/transform_downsampling_representation.R",
    params = params,
    dependencies = list(data),
    duration = duration, memory = "80GB", n_cpus = n_cpus
  )
}

calculate_sirna_quantitative_recovery <- function(dataset, representations, dataset_name, seed, pca_dim, transformations, alphas, regression_knn){
  params <- list(
    input_dataset_id = dataset$result_id,
    representation_result_ids = map_chr(representations, "result_id"),
    dataset = dataset_name,
    seed = seed,
    pca_dim = pca_dim,
    transformations = transformations,
    alphas = alphas,
    regression_knn = regression_knn
  )
  wrap_script(
    "src/downsampling_benchmark/calculate_sirna_quantitative_recovery.R",
    params = params,
    dependencies = c(list(dataset), representations),
    duration = "0-06:00:00", memory = "80GB", n_cpus = 1
  )
}

gather_sirna_quantitative_recovery_results <- function(result_jobs){
  params <- list(sirna_quantitative_recovery_results = map_chr(result_jobs, "result_id"))
  wrap_script(
    "src/downsampling_benchmark/gather_sirna_quantitative_recovery.R",
    params = params,
    dependencies = result_jobs,
    duration = "0-01:00:00", memory = "40GB", n_cpus = 1
  )
}

sirna_quantitative_recovery_jobs <- list()
for(dataset_name in di4$input_data$dataset){
  for(seed in di4$input_data$seed){
    dataset <- get_data_for_downsampling_benchmark(dataset_name, seed)
    for(pca_dim in di4$representation_construction$pca){
      representations <- list()
      for(trans in di4$representation_construction$transformations){
        stopifnot(length(trans$alpha) == 1)
        reduced_rep <- transform_deeply_sequenced_representation(
          trans$name, dataset, data_mode = "reduced", pca_dim = pca_dim, alpha = trans$alpha
        )
        representations <- append(representations, list(reduced_rep))
      }
      res <- calculate_sirna_quantitative_recovery(
        dataset, representations,
        dataset_name = dataset_name, seed = seed, pca_dim = pca_dim,
        transformations = vapply(di4$representation_construction$transformations, function(e) e$name, character(1L)),
        alphas = vapply(di4$representation_construction$transformations, function(e) e$alpha, character(1L)),
        regression_knn = di4$representation_construction$regression_knn[[1]]
      )
      sirna_quantitative_recovery_jobs <- append(sirna_quantitative_recovery_jobs, list(res))
    }
  }
  message("Added all jobs for ", dataset_name, " (siRNA quantitative recovery)")
}

sirna_quantitative_recovery <- gather_sirna_quantitative_recovery_results(sirna_quantitative_recovery_jobs)
table(sapply(sirna_quantitative_recovery$dependencies, job_status))
saveRDS(sirna_quantitative_recovery, file.path("output/job_storage", "sirna_quantitative_recovery_job.RDS"))
sirna_quantitative_recovery <- run_job(sirna_quantitative_recovery, "normal")

if(job_status(sirna_quantitative_recovery) == "done"){
  dir.create("output/benchmark_results", showWarnings = FALSE, recursive = TRUE)
  file.copy(
    file.path(.OUTPUT_FOLDER, "results", sirna_quantitative_recovery$result_id),
    "output/benchmark_results/sirna_quantitative_recovery_results.tsv",
    overwrite = TRUE
  )
  message("Saved output/benchmark_results/sirna_quantitative_recovery_results.tsv")
}else{
  message("Jobs submitted. Re-run this script later to collect the final TSV once dependencies are done.")
}
